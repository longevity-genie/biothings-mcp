#!/usr/bin/env python3
"""Biothings MCP Server - Database query interface for biological data."""

import os
from pathlib import Path
from typing import List, Dict, Any, Optional
import sys

import typer
from fastmcp import FastMCP
from pydantic import BaseModel, Field
from eliot import start_action

from biothings_mcp.biothings_api import GeneTools, VariantTools, ChemTools, TaxonTools
from biothings_mcp.download_api import DownloadTools

# Configuration
DEFAULT_HOST = os.getenv("MCP_HOST", "0.0.0.0")
DEFAULT_PORT = int(os.getenv("MCP_PORT", "3001"))
DEFAULT_TRANSPORT = os.getenv("MCP_TRANSPORT", "streamable-http")
DEFAULT_OUTPUT_DIR = os.getenv("MCP_OUTPUT_DIR", "biothings_output")

class BiothingsMCP(FastMCP):
    """Biothings MCP Server with biological data tools that can be inherited and extended."""
    
    def __init__(
        self, 
        name: str = "Biothings MCP Server",
        prefix: str = "biothings_",
        output_dir: Optional[str] = None,
        **kwargs
    ):
        """Initialize the Biothings tools with FastMCP functionality."""
        # Initialize FastMCP with the provided name and any additional kwargs
        super().__init__(name=name, **kwargs)
        
        self.prefix = prefix
        self.output_dir = output_dir or DEFAULT_OUTPUT_DIR
        self.download_tools = None  # Initialize as None, will be created if needed
        # Register our core tools
        self._register_biothings_tools()
    
    def _register_biothings_tools(self, include_download_tools: bool = False):
        """Register Biothings-specific tools."""
        # Initialize core tool handlers
        self.gene_tools = GeneTools(self, self.prefix)
        self.variant_tools = VariantTools(self, self.prefix)
        self.chem_tools = ChemTools(self, self.prefix)
        self.taxon_tools = TaxonTools(self, self.prefix)
        
        # Register tools from each core handler
        self.gene_tools.register_tools()
        self.variant_tools.register_tools()
        self.chem_tools.register_tools()
        self.taxon_tools.register_tools()
        
        # Conditionally register download tools
        if include_download_tools:
            self.download_tools = DownloadTools(self, self.prefix, self.output_dir)
            self.download_tools.register_tools()
    
    def register_stdio_only(self):
        """Register tools that are only available for stdio transport."""
        if self.download_tools is None:
            self.download_tools = DownloadTools(self, self.prefix, self.output_dir)
            self.download_tools.register_tools()
    
    def run(self, transport: str = "streamable-http", **kwargs):
        """Run the MCP server with conditional download tools registration."""
        # Register stdio-only tools for stdio transport
        if transport == "stdio":
            self.register_stdio_only()
        
        # Call parent run method
        super().run(transport=transport, **kwargs)

# Initialize the Biothings MCP server (will be created per command)
mcp = None

# Create typer app
app = typer.Typer(help="Biothings MCP Server - Database query interface for biological data")

@app.command("run")
def cli_app(
    host: str = typer.Option(DEFAULT_HOST, "--host", help="Host to bind to"),
    port: int = typer.Option(DEFAULT_PORT, "--port", help="Port to bind to"),
    transport: str = typer.Option("streamable-http", "--transport", help="Transport type"),
    output_dir: str = typer.Option(DEFAULT_OUTPUT_DIR, "--output-dir", help="Output directory for local files")
) -> None:
    """Run the MCP server with specified transport."""
    mcp = BiothingsMCP(output_dir=output_dir)
    mcp.run(transport=transport, host=host, port=port)

@app.command("stdio")
def cli_app_stdio(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
    output_dir: str = typer.Option(DEFAULT_OUTPUT_DIR, "--output-dir", help="Output directory for local files")
) -> None:
    """Run the MCP server with stdio transport."""
    mcp = BiothingsMCP(output_dir=output_dir)
    mcp.run(transport="stdio")

@app.command("sse")
def cli_app_sse(
    host: str = typer.Option(DEFAULT_HOST, "--host", help="Host to bind to"),
    port: int = typer.Option(DEFAULT_PORT, "--port", help="Port to bind to"),
    output_dir: str = typer.Option(DEFAULT_OUTPUT_DIR, "--output-dir", help="Output directory for local files")
) -> None:
    """Run the MCP server with SSE transport."""
    mcp = BiothingsMCP(output_dir=output_dir)
    mcp.run(transport="sse", host=host, port=port)

# Standalone CLI functions for direct script access
def cli_app_run() -> None:
    """Standalone function for biothings-mcp-run script."""
    mcp = BiothingsMCP(output_dir=DEFAULT_OUTPUT_DIR)
    mcp.run(transport="streamable-http", host=DEFAULT_HOST, port=DEFAULT_PORT)

def cli_app_stdio() -> None:
    """Standalone function for biothings-mcp-stdio script."""
    mcp = BiothingsMCP(output_dir=DEFAULT_OUTPUT_DIR)
    mcp.run(transport="stdio")

def cli_app_sse() -> None:
    """Standalone function for biothings-mcp-sse script."""
    mcp = BiothingsMCP(output_dir=DEFAULT_OUTPUT_DIR)
    mcp.run(transport="sse", host=DEFAULT_HOST, port=DEFAULT_PORT)

if __name__ == "__main__":
    from pycomfort.logging import to_nice_stdout, to_nice_file
    
    to_nice_stdout()
    # Determine project root and logs directory
    project_root = Path(__file__).resolve().parents[2]
    log_dir = project_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists

    # Define log file paths
    json_log_path = log_dir / "mcp_server.log.json"
    rendered_log_path = log_dir / "mcp_server.log"
    
    # Configure file logging
    to_nice_file(output_file=json_log_path, rendered_file=rendered_log_path)
    app()