# biothings-mcp
[![Tests](https://github.com/longevity-genie/biothings-mcp/actions/workflows/tests.yml/badge.svg)](https://github.com/longevity-genie/biothings-mcp/actions/workflows/tests.yml)
[![PyPI version](https://badge.fury.io/py/biothings-mcp.svg)](https://badge.fury.io/py/biothings-mcp)

MCP (Model Context Protocol) server for biothings

This server implements the Model Context Protocol (MCP) for BioThings, providing a standardized interface for accessing and manipulating biomedical data. MCP enables AI assistants and agents to access specialized biomedical knowledge through structured interfaces to authoritative data sources.

## About MCP (Model Context Protocol)

MCP is a protocol that bridges the gap between AI systems and specialized domain knowledge. It enables:

- **Structured Access**: Direct connection to authoritative biomedical data sources
- **Natural Language Queries**: Simplified interaction with specialized databases
- **Type Safety**: Strong typing and validation through biothings-typed-client
- **AI Integration**: Seamless integration with AI assistants and agents

## Features

- **Type-Safe API**: Built on top of biothings-typed-client for type-safe access to BioThings services
- **Standardized Protocol**: Implements the MCP specification for consistent data access
- **Multiple Data Types**: Supports access to various biomedical data types including:
  - Genes
  - Variants
  - Chemicals
  - Taxons
- **Synchronous & Asynchronous**: Supports both synchronous and asynchronous operations
- **Validation**: Built-in data validation and type checking
- **AI Integration**: Designed for seamless integration with AI systems and agents

## Available Clients

The server provides access to the following typed clients from biothings-typed-client:

- `VariantClient` / `VariantClientAsync`: For accessing variant data
- `GeneClient` / `GeneClientAsync`: For accessing gene data
- `ChemClient` / `ChemClientAsync`: For accessing chemical compound data
- `TaxonClient` / `TaxonClientAsync`: For accessing taxonomic information

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
# Start the MCP server
uv run server
```

### Integration with AI Systems

Configure your AI system to use the MCP server in one of two ways:

1. Direct SSE Connection:
```json
{
  "mcpServers": {
    "biothings-mcp": {
      "url": "http://localhost:3001/mcp"
    }
  }
}
```

2. Using mcp-remote (recommended for OAuth support):
```json
{
  "mcpServers": {
    "biothings-mcp": {
      "command": "npx",
      "args": [
        "mcp-remote",
        "http://localhost:3001/mcp",
        "8080"  // Optional port number for OAuth support
      ]
    }
  }
}
```

## Testing & Verification

Run tests for the API endpoint:
```bash
uv run pytest -vvv -s
```

Test your MCP setup with the MCP Inspector:

```bash
npx @modelcontextprotocol/inspector --config mcp-config.json --server biothings-mcp
```

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
