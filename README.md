# biothings-mcp
[![Tests](https://github.com/longevity-genie/biothings-mcp/actions/workflows/tests.yml/badge.svg)](https://github.com/longevity-genie/biothings-mcp/actions/workflows/tests.yml)
[![PyPI version](https://badge.fury.io/py/biothings-mcp.svg)](https://badge.fury.io/py/biothings-mcp)

MCP (Model Context Protocol) server for Biothings.io

This server implements the Model Context Protocol (MCP) for BioThings, providing a standardized interface for accessing and manipulating biomedical data. MCP enables AI assistants and agents to access specialized biomedical knowledge through structured interfaces to authoritative data sources. Supported BioThings data sources include:

- [mygene.info](https://mygene.info) — Gene annotation and query service
- [myvariant.info](https://myvariant.info) — Variant annotation and query service
- [mychem.info](https://mychem.info) — Chemical compound annotation and query service

If you want to understand more what is Model Context Protocol and how to use it more efficiently you can take [DeepLearning AI Course](https://www.deeplearning.ai/short-courses/mcp-build-rich-context-ai-apps-with-anthropic/) or just search for MCP videos on YouTube.


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
uvx --version
```

uvx is a very nice tool that can run a python package installing it if needed.

### Running with uvx

You can run the biothings-mcp server directly using uvx without cloning the repository:

#### STDIO Mode (for MCP clients that require stdio, can be useful when you want to save files)
```bash
# Run the server in STDIO mode (default mode)
uvx biothings-mcp

# Or explicitly specify stdio mode
uvx --from biothings-mcp stdio
```

#### HTTP Mode (Web Server)
```bash
# Run the server in streamable HTTP mode on default ( 3001 ) port
uvx biothings-mcp

# Run on a specific port
uvx biothings-mcp server --port 8000
```

#### SSE Mode (Server-Sent Events)
```bash
# Run the server in SSE mode
uvx biothings-mcp sse
```

The HTTP mode will start a web server that you can access at `http://localhost:3001/mcp` (with documentation at `http://localhost:3001/docs`). The STDIO mode is designed for MCP clients that communicate via standard input/output, while SSE mode uses Server-Sent Events for real-time communication.


## Configuring your (Anthropic Claude Desktop, Cursor, Windsurf, etc.)

We provide stdio configuration using the proxy (might need npx to run):
* mcp-config-remote.json - for remote configuration
* mcp-config-stdio.json - stdio configuration for localhost for MCP clients which do not support 


### Inspecting Biothings MCP server

If you want to inspect the methods provided by the MCP use npx (you may need to install nodejs and npm)

Test your MCP setup with the MCP Inspector.

If you want to inspect local streamable-http server you use:

```bash
npx @modelcontextprotocol/inspector --config mcp-config.json --server biothings-mcp
```

Add -remote suffix for the remote server.


If you want to inspect stdio local server you should use

```bash
npx @modelcontextprotocol/inspector --config mcp-config-stdio.json --server biothings-mcp
```

You can also run inspector manually and put server parameters in the interface:
```
npx @modelcontextprotocol/inspector
```

After that you can explore its methods with MCP Inspector at http://127.0.0.1:6274


## Repository setup

```bash
# Clone the repository
git clone git@github.com:longevity-genie/biothings-mcp.git
cd biothings-mcp
uv sync
```

### Running the MCP Server

If you already cloned the repo you can run the server with uv

```bash
# Start the MCP server locally
uv run server
```


### Integration with AI Systems

To integrate this server with your MCP-compatible AI client, you can use one of the preconfigured JSON files provided in this repository:

*   **For connecting to a locally running server:** Use `mcp-config.json`. Ensure the server is running first, either via `uv run server` (see [Running the MCP Server](#running-the-mcp-server)) or `docker-compose up` (see [Docker Deployment](#docker-deployment)).
*   **For connecting to the publicly hosted server:** Use `mcp-config-remote.json`. This connects to `https://biothings.longevity-genie.info/mcp` and doesn't require you to run anything locally.

Simply point your AI client (like Cursor, Windserve, ClaudeDesktop, VS Code with Copilot, or [others](https://github.com/punkpeye/awesome-mcp-clients)) to use the appropriate configuration file.

Here's an example of how the tools might appear in an MCP client like Cursor after configuration:

![Cursor Usage Example](images/cursor_usage_example.jpg)

## KNOWN ISSUES

The library is beta-quality. The major problem right now is that LLM-s are often stupid and do not know how to put valid gene and gene variant symbols. We plan to mitigrate it by extending comments and providing additional method for entity resolution.

## Testing & Verification

Run tests for the API endpoint:
```bash
uv run pytest -vvv -s
```

You can use MCP inspector with locally build MCP server same way as with uvx

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
