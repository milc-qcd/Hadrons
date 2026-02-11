#!/usr/bin/env bash

USE_CUSTOM=false
for arg in "$@"; do
    if [ "$arg" = "--custom" ]; then
        USE_CUSTOM=true
        break
    fi
done

if [ "$USE_CUSTOM" = true ]; then
    rm -f modules.inc
    ln -s modules_custom.inc modules.inc
    MODULES_FILE="modules_custom.inc"
else
    MODULES_FILE="modules.inc"

    echo 'modules_cpp =\' > modules.inc
    find Modules -name '*.cpp' -type f -print | LC_ALL=C sort | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
    echo '' >> modules.inc
    echo 'modules_hpp =\' >> modules.inc
    find Modules -name '*.hpp' -type f -print | LC_ALL=C sort | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
    echo '' >> modules.inc
fi

# Generate Modules.hpp from modules_hpp in MODULES_FILE
rm -f Modules.hpp
# Extract hpp files from modules_hpp variable, removing leading spaces and trailing backslashes
sed -n '/^modules_hpp/,/^$/p' "$MODULES_FILE" | grep '\.hpp' | sed 's/^[[:space:]]*//;s/[[:space:]]*\\$//' | while read -r f; do
	if [ -n "$f" ]; then
		echo "#include <Hadrons/${f}>" >> Modules.hpp
	fi
done
