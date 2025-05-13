FROM ghcr.io/astral-sh/uv:debian-slim

# Set the working directory in the container
WORKDIR /app

COPY . .

RUN uv sync

RUN uv run pytest -vvv

# Expose the port the app runs on (adjust if your server runs on a different port)
EXPOSE 8000

# Define the command to run the application
CMD ["uv", "run", "server"] 