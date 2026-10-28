#!/usr/bin/env bash

set -euo pipefail

if [ "$(uname -s)" != "Darwin" ]; then
    echo "binary-mac must be run on macOS." >&2
    exit 1
fi

command -v otool >/dev/null || { echo "otool is required (install Xcode Command Line Tools)." >&2; exit 1; }
command -v install_name_tool >/dev/null || { echo "install_name_tool is required (install Xcode Command Line Tools)." >&2; exit 1; }
command -v codesign >/dev/null || { echo "codesign is required (install Xcode Command Line Tools)." >&2; exit 1; }

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
package_name="flick-macos-$arch"
dist_dir="$project_dir/dist"
stage="$dist_dir/$package_name"
archive="$dist_dir/$package_name.tar.gz"
mapping=$(mktemp "${TMPDIR:-/tmp}/flick-libraries.XXXXXX")
trap 'rm -f "$mapping"' EXIT

rm -rf "$stage"
mkdir -p "$stage/bin" "$stage/lib" "$stage/resources/external" "$dist_dir"
cp "$executable" "$stage/bin/flick-bin"

# Retain the complete project data set.  Source and build-only files are
# excluded, except for the explicitly requested Eigen and AccuRT trees.
for item in astronomy coating component environment geometry material mie model numeric polarization radiator transporter main/commands Example; do
    if [ -e "$project_dir/$item" ]; then
        mkdir -p "$stage/resources/$(dirname "$item")"
        rsync -a --exclude='*.o' --exclude='*.d' --exclude='test_all' --exclude='flick' \
            "$project_dir/$item" "$stage/resources/$(dirname "$item")/"
    fi
done
rsync -a --exclude='.git' "$project_dir/external/eigen" "$stage/resources/external/"
rsync -a --exclude='.git' "$project_dir/external/accurt-kernel" "$stage/resources/external/"
cp "$project_dir/README.md" "$stage/resources/"

# Copy every non-Apple dynamic dependency. Repeat because libraries can depend
# on further Homebrew/GCC libraries of their own.
added_libraries=1
while [ "$added_libraries" -eq 1 ]; do
    added_libraries=0
    while IFS= read -r -d '' scan_file; do
        while IFS= read -r dependency; do
            case "$dependency" in
                /System/*|/usr/lib/*|@executable_path/*|@loader_path/*) continue ;;
                @rpath/*)
                    dependency_name=${dependency#@rpath/}
                    scan_name=$(basename "$scan_file")
                    scan_origin=$(awk -F '|' -v name="$scan_name" '$2 == name { print $1; exit }' "$mapping")
                    candidate=$(dirname "$scan_origin")/$dependency_name
                    if [ -n "$scan_origin" ] && [ -f "$candidate" ]; then
                        dependency=$candidate
                    elif [ -f "$stage/lib/$dependency_name" ]; then
                        continue
                    else
                        echo "Cannot locate unresolved dependency '$dependency' required by $scan_file" >&2
                        exit 1
                    fi ;;
            esac
            [ -f "$dependency" ] || { echo "Missing dynamic library: $dependency" >&2; exit 1; }
            library_name=$(basename "$dependency")
            destination="$stage/lib/$library_name"
            if [ ! -e "$destination" ]; then
                cp "$dependency" "$destination"
                chmod u+w "$destination"
                printf '%s|%s\n' "$dependency" "$library_name" >> "$mapping"
                added_libraries=1
            else
                recorded=$(awk -F '|' -v name="$library_name" '$2 == name { print $1; exit }' "$mapping")
                if [ -n "$recorded" ] && [ "$recorded" != "$dependency" ]; then
                    echo "Two dependencies have the same filename: $recorded and $dependency" >&2
                    exit 1
                fi
            fi
        done < <(otool -L "$scan_file" | tail -n +2 | sed -E 's/^[[:space:]]*([^[:space:]]+).*/\1/')
    done < <(find "$stage/bin" "$stage/lib" -type f -print0)
done

# Point the executable and bundled libraries at the package-local lib folder.
for binary in "$stage/bin/flick-bin" "$stage"/lib/*.dylib; do
    [ -e "$binary" ] || continue
    while IFS='|' read -r original name; do
        [ -n "$original" ] || continue
        install_name_tool -change "$original" "@loader_path/../lib/$name" "$binary" 2>/dev/null || true
        install_name_tool -change "@rpath/$name" "@loader_path/../lib/$name" "$binary" 2>/dev/null || true
    done < "$mapping"
done
for library in "$stage"/lib/*.dylib; do
    [ -e "$library" ] || continue
    install_name_tool -id "@rpath/$(basename "$library")" "$library"
    codesign --force --sign - "$library" >/dev/null
done
codesign --force --sign - "$stage/bin/flick-bin" >/dev/null

cat > "$stage/flick" <<'LAUNCHER'
#!/usr/bin/env bash
set -euo pipefail

# Follow an installed symlink back to the self-contained Flick package.
launcher=$0
while [ -L "$launcher" ]; do
    launcher_dir=$(cd -P -- "$(dirname -- "$launcher")" && pwd -P)
    launcher=$(readlink "$launcher")
    case "$launcher" in
        /*) ;;
        *) launcher="$launcher_dir/$launcher" ;;
    esac
done

package_dir=$(cd -P -- "$(dirname -- "$launcher")" && pwd -P)
export FLICK_PATH="$package_dir/resources"
export DYLD_LIBRARY_PATH="$package_dir/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
exec "$package_dir/bin/flick-bin" "$@"
LAUNCHER
chmod +x "$stage/flick"

cat > "$stage/install.sh" <<'INSTALLER'
#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./install.sh [--prefix DIRECTORY] [--force]

Install Flick for the current user (default prefix: $HOME/.local).
Use --force to replace an existing installation at the selected prefix.
EOF
}

prefix=${HOME:?HOME is not set}/.local
force=0
while [ "$#" -gt 0 ]; do
    case "$1" in
        --prefix)
            [ "$#" -ge 2 ] || { echo "--prefix requires a directory" >&2; exit 2; }
            prefix=$2
            shift 2
            ;;
        --force)
            force=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$prefix" in
    /*) ;;
    *) prefix=$(cd -P -- "$(dirname -- "$prefix")" && pwd -P)/$(basename -- "$prefix") ;;
esac

source_dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
install_dir="$prefix/share/flick"
bin_dir="$prefix/bin"
launcher="$bin_dir/flick"

if [ "$source_dir" = "$install_dir" ]; then
    echo "Flick is already installed in $install_dir"
    exit 0
fi
if [ -e "$install_dir" ] && [ "$force" -ne 1 ]; then
    echo "Flick is already installed in $install_dir" >&2
    echo "Run with --force to replace it." >&2
    exit 1
fi
if { [ -e "$launcher" ] || [ -L "$launcher" ]; } && \
   { [ ! -L "$launcher" ] || [ "$(readlink "$launcher")" != "$install_dir/flick" ]; }; then
    echo "Refusing to replace $launcher because it does not belong to this installation." >&2
    exit 1
fi

mkdir -p "$prefix/share" "$bin_dir"
temporary="$prefix/share/.flick-install.$$"
backup="$prefix/share/.flick-backup.$$"
cleanup() {
    rm -rf "$temporary"
    if [ -e "$backup" ]; then
        if [ ! -e "$install_dir" ]; then
            mv "$backup" "$install_dir"
        else
            rm -rf "$backup"
        fi
    fi
}
trap cleanup EXIT
cp -R "$source_dir" "$temporary"

if [ -e "$install_dir" ]; then
    mv "$install_dir" "$backup"
fi
mv "$temporary" "$install_dir"
ln -sfn "$install_dir/flick" "$launcher"
rm -rf "$backup"
trap - EXIT

echo "Installed Flick in $install_dir"
echo "Created $launcher"

case ":${PATH:-}:" in
    *:"$bin_dir":*) ;;
    *)
        echo
        echo "$bin_dir is not currently in PATH. Add this line to your shell profile:"
        echo "  export PATH=\"$bin_dir:\$PATH\""
        ;;
esac
INSTALLER
chmod +x "$stage/install.sh"

cat > "$stage/uninstall.sh" <<'UNINSTALLER'
#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: ./uninstall.sh [--prefix DIRECTORY]

Remove Flick from the selected prefix (default: $HOME/.local).
EOF
}

prefix=${HOME:?HOME is not set}/.local
while [ "$#" -gt 0 ]; do
    case "$1" in
        --prefix)
            [ "$#" -ge 2 ] || { echo "--prefix requires a directory" >&2; exit 2; }
            prefix=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$prefix" in
    /*) ;;
    *) prefix=$(cd -P -- "$(dirname -- "$prefix")" && pwd -P)/$(basename -- "$prefix") ;;
esac

install_dir="$prefix/share/flick"
launcher="$prefix/bin/flick"
expected="$install_dir/flick"

if [ -L "$launcher" ] && [ "$(readlink "$launcher")" = "$expected" ]; then
    rm "$launcher"
elif [ -e "$launcher" ] || [ -L "$launcher" ]; then
    echo "Refusing to remove $launcher because it does not belong to this installation." >&2
    exit 1
fi

if [ -d "$install_dir" ]; then
    rm -rf "$install_dir"
    echo "Removed $install_dir"
else
    echo "Flick is not installed in $install_dir"
fi
UNINSTALLER
chmod +x "$stage/uninstall.sh"

cat > "$stage/README.txt" <<EOF
Flick for macOS ($arch)

Run from Terminal:

    ./flick help

Install for the current user so that it can be run from any directory:

    ./install.sh

If ~/.local/bin is not already in PATH, the installer prints the line to add
to your shell profile. To select another location, use for example:

    ./install.sh --prefix /usr/local

Remove that installation with:

    ./uninstall.sh

This package is self-contained and does not require compilation. It is built
for $arch Macs and may require a macOS version at least as new as the build Mac.
The AccuRT kernel is distributed under its license in:

    resources/external/accurt-kernel/AccuRT-license-NC-v2.pdf
EOF

for script in "$stage/flick" "$stage/install.sh" "$stage/uninstall.sh"; do
    bash -n "$script"
done
"$stage/flick" help >/dev/null
tar -C "$dist_dir" -czf "$archive" "$package_name"
echo "Created $archive"
