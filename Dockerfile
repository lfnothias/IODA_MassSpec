# Use an appropriate base image with Python preinstalled, e.g., python:3.10
FROM python:3.10

# Copy your wheel files to the container
COPY pyopenms_wheels/ /pyopenms_wheels/

# Run the installation script
RUN set -ex; \
    PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')"); \
    echo "Detected Python version: ${PYTHON_VERSION}"; \
    PLATFORM_SYSTEM=$(uname); \
    echo "Detected platform: ${PLATFORM_SYSTEM}"; \
    if [ "$PLATFORM_SYSTEM" == "Windows" ]; then \
        WHL_FILE="/pyopenms_wheels/pyopenms-3.0.0.dev20230316-cp${PYTHON_VERSION}-cp${PYTHON_VERSION}-win_amd64.whl"; \
    elif [ "$PLATFORM_SYSTEM" == "Darwin" ]; then \
        WHL_FILE="/pyopenms_wheels/pyopenms-3.0.0.dev20230316-cp${PYTHON_VERSION}-cp${PYTHON_VERSION}-macosx_10_9


# Set working directory
WORKDIR /app

# Copy requirements.txt and other necessary files
COPY requirements.txt /app/requirements.txt

# Install requirements
RUN pip install --no-cache-dir -r requirements.txt
