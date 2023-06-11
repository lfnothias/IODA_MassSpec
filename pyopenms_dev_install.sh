#!/bin/bash

# Install other required packages if needed

PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")
echo "Detected Python version: ${PYTHON_VERSION}"

PLATFORM_SYSTEM=$(uname)
echo "Detected platform: ${PLATFORM_SYSTEM}"

if [ "$PLATFORM_SYSTEM" == "Windows" ]; then
    WHL_URL="https://github.com/OpenMS/OpenMS/suites/13420941794/artifacts/735505003"
elif [ "$PLATFORM_SYSTEM" == "Darwin" ]; then
    WHL_URL="https://github.com/OpenMS/OpenMS/suites/13420941794/artifacts/735505002"
else
    WHL_URL="https://github.com/OpenMS/OpenMS/suites/13420941794/artifacts/735505000"
fi

echo "Selected wheel URL: ${WHL_URL}"
WHEEL_FILE=$(basename "${WHL_URL}")
curl -OL $WHL_URL
pip install $WHEEL_FILE
