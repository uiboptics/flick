#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Flick environment setup script
# Detects shell, sets FLICK_PATH and compiler, and updates shell profile safely.
# -----------------------------------------------------------------------------

set -euo pipefail

FLICK_FOLDER=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

# Detect shell and select appropriate profile + compiler
if [ -n "$($SHELL -c 'echo $ZSH_VERSION')" ]; then
    PROFILE="$HOME/.zprofile"
    COMPILER='clang++ -std=c++20'
elif [ -n "$($SHELL -c 'echo $BASH_VERSION')" ]; then
    PROFILE="$HOME/.bashrc"
    COMPILER='g++ -std=c++20'
else
    echo "⚠️  Unknown shell. Please set environmental variables manualla"
    echo "     based on information in update_shell.sh."
    exit 1
fi

# Check if FLICK_PATH already exists in the profile
if grep -q 'FLICK_PATH' "$PROFILE"; then
    echo ""
    echo "ℹ️  Flick environment variables already found in $PROFILE."
    echo "    No changes made."
    echo ""
else
    {
        echo ""
        echo "# --- Flick environment variables ---"
	echo "export FLICK_PATH=$FLICK_FOLDER"
	echo "export FLICK_COMPILER=\"$COMPILER\""
	echo "export EIGEN_PATH=\$FLICK_PATH/external/eigen"
        echo "export CPLUS_INCLUDE_PATH=\$FLICK_PATH:\$EIGEN_PATH:\$CPLUS_INCLUDE_PATH"
        echo "export PATH=\$FLICK_PATH/main:\$PATH"
        echo "# --- End Flick environment variables ---"
    } >> "$PROFILE"

    echo ""
    echo "✅ Flick environment variables added to $PROFILE"
    echo ""
    echo "Now run:"
    echo
    echo "  source $PROFILE"
    echo
    echo "then:"
    echo
    echo "  make"
    echo
fi
