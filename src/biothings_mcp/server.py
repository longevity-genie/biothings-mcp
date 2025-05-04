import os
from pathlib import Path
from typing import Optional, List, Dict, Any, Union, Type, ClassVar

from dotenv import load_dotenv
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi_mcp import FastApiMCP

from biothings_mcp.biothings_api import BiothingsRestAPI
from pycomfort.logging import to_nice_stdout, to_nice_file
# from biothings_mcp.logging import configure_logging, LoggedTask, log_info

# Load environment variables
load_dotenv()

# Configuration
DEFAULT_HOST = os.getenv("MCP_HOST", "0.0.0.0")
DEFAULT_PORT = int(os.getenv("MCP_PORT", "3001"))

def create_app() -> FastAPI:
    """Create and configure the FastAPI application"""
    # Configure logging
    # configure_logging()
    # log_info("Starting Biothings MCP server")
    
    # with LoggedTask("create_app") as task:
    app = BiothingsRestAPI()
        
    # Configure CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
        
    # Mount MCP
    mcp = FastApiMCP(app)
    mcp.mount()
        
    # log_info("Biothings MCP server configured successfully")
    return app

def run_server(host: str = DEFAULT_HOST, port: int = DEFAULT_PORT):
    """Entry point for running the server"""
    import uvicorn
    # with LoggedTask("run_server", host=host, port=port) as task:
    app = create_app()
    # log_info(f"Starting server on {host}:{port}")
    uvicorn.run(app, host=host, port=port)

if __name__ == "__main__":
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
    run_server()