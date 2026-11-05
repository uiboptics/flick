#!/usr/bin/env bash

set -euo pipefail

if [ "$(uname -s)" != "Linux" ]; then
    echo "binary-linux must be run on Linux." >&2
    exit 1
fi

command -v ldd >/dev/null || { echo "ldd is required." >&2; exit 1; }
command -v rsync >/dev/null || { echo "rsync is required." >&2; exit 1; }

project_dir=$(cd -P -- "$(dirname -- "$0")/.." && pwd -P)
executable="$project_dir/main/flick"

if [ ! -x "$executable" ]; then
    echo "The Flick executable was not built: $executable" >&2
    exit 1
fi
if [ ! -d "$project_dir/external/eigen/Eigen" ]; then
    echo "Eigen is missing. Run 'make install-eigen'." >&2
    exit 1
fi
if [ ! -f "$project_dir/external/accurt-kernel/kernel.f95" ]; then
    echo "The AccuRT kernel is missing. Run 'make install-accurt-kernel'." >&2
    exit 1
fi

arch=$(uname -m)
package_name="flick-linux-$arch"
dist_dir="$project_dir/dist"
stage="$dist_dir/$package_name"
archive="$dist_dir/$package_name.tar.gz"

rm -rf "$stage"
mkdir -p "$stage/bin" "$stage/lib" "$stage/python-package/python_script" \
    "$stage/resources/external" "$dist_dir"
cp "$executable" "$stage/bin/flick-bin"

cp "$project_dir/pyproject.toml" "$project_dir/README.md" \
    "$stage/python-package/"
cp "$project_dir/python_script/flick.py" \
    "$project_dir/python_script/snow_optics.py" \
    "$stage/python-package/python_script/"

for item in astronomy coating component environment geometry material mie model numeric polarization radiator transporter main/commands Example; do
    if [ -e "$project_dir/$item" ]; then
        mkdir -p "$stage/resources/$(dirname "$item")"
        rsync -a --exclude='*.o' --exclude='*.d' --exclude='test_all' --exclude='flick' \
            "$project_dir/$item" "$stage/resources/$(dirname "$item")/"
    fi
done
# Eigen and the AccuRT Fortran source are build-time inputs and must not be
# distributed. Keep the AccuRT license alongside the compiled-in kernel.
mkdir -p "$stage/resources/external/accurt-kernel"
cp "$project_dir/external/accurt-kernel/AccuRT-license-NC-v2.pdf" \
    "$stage/resources/external/accurt-kernel/"
cp "$project_dir/README.md" "$stage/resources/"

# Fail closed if either build-time source tree is accidentally staged later.
[ ! -e "$stage/resources/external/eigen" ] || {
    echo "Refusing to package the Eigen source tree." >&2
    exit 1
}
[ ! -e "$stage/resources/external/accurt-kernel/kernel.f95" ] || {
    echo "Refusing to package the AccuRT kernel source." >&2
    exit 1
}

# Bundle compiler and third-party libraries, but leave glibc and the dynamic
# loader to the target system. LD_LIBRARY_PATH in the launcher also lets these
# copied libraries find one another without modifying the ELF files.
copy_dependencies() {
    local scan_file=$1 dependency name dependencies
    dependencies=$(ldd "$scan_file")
    if grep -q '=> not found' <<<"$dependencies"; then
        grep '=> not found' <<<"$dependencies" >&2
        return 1
    fi
    while IFS= read -r dependency; do
        [ -f "$dependency" ] || continue
        name=$(basename "$dependency")
        case "$name" in
            ld-linux*|libc.so.*|libdl.so.*|libm.so.*|libpthread.so.*|libresolv.so.*|librt.so.*|libutil.so.*) continue ;;
        esac
        if [ ! -e "$stage/lib/$name" ]; then
            cp -L "$dependency" "$stage/lib/$name"
            copy_dependencies "$stage/lib/$name"
        fi
    done < <(awk '
        /=> \/[^ ]+/ { print $3 }
        /^[[:space:]]*\// { print $1 }
    ' <<<"$dependencies")
}
copy_dependencies "$stage/bin/flick-bin"

cat > "$stage/flick" <<'LAUNCHER'
#!/usr/bin/env bash
set -euo pipefail
launcher=$0
while [ -L "$launcher" ]; do
    launcher_dir=$(cd -P -- "$(dirname -- "$launcher")" && pwd -P)
    launcher=$(readlink "$launcher")
    case "$launcher" in /*) ;; *) launcher="$launcher_dir/$launcher" ;; esac
done
package_dir=$(cd -P -- "$(dirname -- "$launcher")" && pwd -P)
export FLICK_PATH="$package_dir/resources"
export LD_LIBRARY_PATH="$package_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$package_dir/bin/flick-bin" "$@"
LAUNCHER
chmod +x "$stage/flick"

cat > "$stage/install.sh" <<'INSTALLER'
#!/usr/bin/env bash
set -euo pipefail
prefix=${HOME:?HOME is not set}/.local
force=0
while [ "$#" -gt 0 ]; do
    case "$1" in
        --prefix) [ "$#" -ge 2 ] || { echo "--prefix requires a directory" >&2; exit 2; }; prefix=$2; shift 2 ;;
        --force) force=1; shift ;;
        -h|--help) echo "Usage: ./install.sh [--prefix DIRECTORY] [--force]"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$prefix" in /*) ;; *) prefix=$(cd -P -- "$(dirname -- "$prefix")" && pwd -P)/$(basename -- "$prefix") ;; esac
source_dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
install_dir="$prefix/share/flick"
launcher="$prefix/bin/flick"
[ ! -e "$install_dir" ] || [ "$force" -eq 1 ] || { echo "Flick is already installed in $install_dir (use --force to replace it)." >&2; exit 1; }
if { [ -e "$launcher" ] || [ -L "$launcher" ]; } && { [ ! -L "$launcher" ] || [ "$(readlink "$launcher")" != "$install_dir/flick" ]; }; then
    echo "Refusing to replace $launcher because it does not belong to this installation." >&2; exit 1
fi
command -v python3 >/dev/null && python3 -m pip --version >/dev/null || { echo "python3 and pip are required." >&2; exit 1; }
mkdir -p "$prefix/share" "$prefix/bin"
temporary="$prefix/share/.flick-install.$$"
trap 'rm -rf "$temporary"' EXIT
cp -R "$source_dir" "$temporary"
rm -rf "$install_dir"
mv "$temporary" "$install_dir"
ln -sfn "$install_dir/flick" "$launcher"
python3 -m pip install --user "$install_dir/python-package"
trap - EXIT
echo "Installed Flick in $install_dir"
echo "Created $launcher"
INSTALLER
chmod +x "$stage/install.sh"

cat > "$stage/uninstall.sh" <<'UNINSTALLER'
#!/usr/bin/env bash
set -euo pipefail
prefix=${HOME:?HOME is not set}/.local
if [ "${1:-}" = "--prefix" ] && [ "$#" -ge 2 ]; then prefix=$2; shift 2; fi
[ "$#" -eq 0 ] || { echo "Usage: ./uninstall.sh [--prefix DIRECTORY]" >&2; exit 2; }
case "$prefix" in /*) ;; *) prefix=$(cd -P -- "$(dirname -- "$prefix")" && pwd -P)/$(basename -- "$prefix") ;; esac
install_dir="$prefix/share/flick"
launcher="$prefix/bin/flick"
if [ -L "$launcher" ] && [ "$(readlink "$launcher")" = "$install_dir/flick" ]; then rm "$launcher"; elif [ -e "$launcher" ] || [ -L "$launcher" ]; then echo "Refusing to remove unrelated $launcher" >&2; exit 1; fi
if command -v python3 >/dev/null && python3 -m pip --version >/dev/null 2>&1; then python3 -m pip uninstall -y uiboptics-flick; fi
rm -rf "$install_dir"
echo "Removed $install_dir"
UNINSTALLER
chmod +x "$stage/uninstall.sh"

cat > "$stage/README.txt" <<EOF
Flick for Linux ($arch)

Run ./flick help, or run ./install.sh to install for the current user.
Use ./install.sh --prefix /usr/local to select another prefix and
./uninstall.sh to remove the installation.

The binary is built for $arch Linux systems using a glibc version at least as
new as the build system. The installer requires python3 and pip for the Flick
Python interface.
EOF

for script in "$stage/flick" "$stage/install.sh" "$stage/uninstall.sh"; do bash -n "$script"; done
"$stage/flick" help >/dev/null
tar -C "$dist_dir" -czf "$archive" "$package_name"
echo "Created $archive"
