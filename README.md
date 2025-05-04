# biothings-mcp
[![Tests](https://github.com/longevity-genie/biothings-mcp/actions/workflows/tests.yml/badge.svg)](https://github.com/longevity-genie/biothings-mcp/actions/workflows/tests.yml)
[![PyPI version](https://badge.fury.io/py/biothings-mcp.svg)](https://badge.fury.io/py/biothings-mcp)

MCP (Model Context Protocol) server for Biothings.io

This server implements the Model Context Protocol (MCP) for BioThings, providing a standardized interface for accessing and manipulating biomedical data. MCP enables AI assistants and agents to access specialized biomedical knowledge through structured interfaces to authoritative data sources. Supported BioThings data sources include:

- [mygene.info](https://mygene.info) — Gene annotation and query service
- [myvariant.info](https://myvariant.info) — Variant annotation and query service
- [mychem.info](https://mychem.info) — Chemical compound annotation and query service

## About MCP (Model Context Protocol)

MCP is a protocol that bridges the gap between AI systems and specialized domain knowledge. It enables:

- **Structured Access**: Direct connection to authoritative biomedical data sources
- **Natural Language Queries**: Simplified interaction with specialized databases
- **Type Safety**: Strong typing and validation through biothings-typed-client
- **AI Integration**: Seamless integration with AI assistants and agents

## Available API Interfaces

This server provides dedicated API interfaces for different BioThings data types, leveraging the `biothings-typed-client` library. These interfaces are implemented using the following mixins:

- **Gene Interface**: `GeneRoutesMixin` (wraps `GeneClientAsync`)
- **Variant Interface**: `VariantsRoutesMixin` (wraps `VariantClientAsync`)
- **Chemical Interface**: `ChemRoutesMixin` (wraps `ChemClientAsync`)
- **Taxon Interface**: `TaxonRoutesMixin` (wraps `TaxonClientAsync`)

## Quick Start

### Installing uv

```bash
# Download and install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Verify installation
uv --version
```

### Setup

```bash
# Clone the repository
git clone git@github.com:longevity-genie/biothings-mcp.git
cd biothings-mcp
uv sync
```

### Running the MCP Server

```bash
# Start the MCP server locally
uv run server
```

### Docker Deployment

The easiest way to run the MCP server is using Docker. The project provides a pre-built Docker image available on GitHub Container Registry.

1. Using Docker Compose (recommended):

```bash
# Clone the repository
git clone git@github.com:longevity-genie/biothings-mcp.git
cd biothings-mcp

# Start the services
docker-compose up
```

This will start:
- The MCP server on port 3001
- The MCP Inspector on port 6277

2. Using Docker directly:

```bash
# Pull the latest image
docker pull ghcr.io/longevity-genie/biothings-mcp:latest

# Run the container
docker run -p 3001:3001 -e MCP_PORT=3001 ghcr.io/longevity-genie/biothings-mcp:latest
```

The MCP server will be available at `http://localhost:3001/mcp` (with docs at http://localhost:3001/docs).

A publicly hosted version of this server is also available at `https://biothings.longevity-genie.info/mcp` (with docs at https://biothings.longevity-genie.info/docs)

### Integration with AI Systems

To integrate this server with your MCP-compatible AI client, you can use one of the preconfigured JSON files provided in this repository:

*   **For connecting to a locally running server:** Use `mcp-config.json`. Ensure the server is running first, either via `uv run server` (see [Running the MCP Server](#running-the-mcp-server)) or `docker-compose up` (see [Docker Deployment](#docker-deployment)).
*   **For connecting to the publicly hosted server:** Use `mcp-config-remote.json`. This connects to `https://biothings.longevity-genie.info/mcp` and doesn't require you to run anything locally.

Simply point your AI client (like Cursor, Windserve, ClaudeDesktop, VS Code with Copilot, or [others](https://github.com/punkpeye/awesome-mcp-clients)) to use the appropriate configuration file.

Here's an example of how the tools might appear in an MCP client like Cursor after configuration:

![Cursor Usage Example](images/cursor_usage_example.jpg)

## KNOWN ISSUES

The library is alpha-quality. The major problem right now is that LLM-s are often stupid and do not know how to put valid gene and gene variant symbols. We plan to mitigrate it by extending comments and providing additional method for entetity resolution.

## Testing & Verification

Run tests for the API endpoint:
```bash
uv run pytest -vvv -s
```

Test your MCP setup with the MCP Inspector:

```bash
npx @modelcontextprotocol/inspector --config mcp-config.json --server biothings-mcp
```
After that you can explore its methods with MCP Inspector at http://127.0.0.1:6274

*Note: Using the MCP Inspector is optional. Most MCP clients (like Cursor, Windsurv, etc.) will automatically display the available tools from this server once configured. However, the Inspector can be useful for detailed testing and exploration.* 

*If you choose to use the Inspector via `npx`, ensure you have Node.js and npm installed. Using [nvm](https://github.com/nvm-sh/nvm) (Node Version Manager) is recommended for managing Node.js versions.*

This opens a web interface where you can explore and test all available tools.

## Documentation

For detailed documentation about the MCP protocol and its implementation, refer to:
- [MCP Protocol Documentation](https://modelcontextprotocol.org)
- [biothings-typed-client Documentation](https://github.com/longevity-genie/biothings-typed-client)
- [FastAPI-MCP Documentation](https://github.com/tadata-org/fastapi_mcp)

## License

This project is licensed under the MIT License.

## Acknowledgments

- [BioThings](https://biothings.io/) for the REST API and original [client library](https://github.com/biothings/biothings_client.py)
- [MCP Protocol](https://modelcontextprotocol.org) for the protocol specification
- [Pydantic](https://pydantic-docs.helpmanual.io/) for the data validation framework
- [FastAPI-MCP](https://github.com/tadata-org/fastapi_mcp) for the MCP server implementation

- This project is part of the [Longevity Genie](https://github.com/longevity-genie) organization, which develops open-source AI assistants and libraries for health, genetics, and longevity research.

We are supported by:

[![HEALES](images/heales.jpg)](https://heales.org/)

*HEALES - Healthy Life Extension Society*

and

[![IBIMA](images/IBIMA.jpg)](https://ibima.med.uni-rostock.de/)

[IBIMA - Institute for Biostatistics and Informatics in Medicine and Ageing Research](https://ibima.med.uni-rostock.de/)
