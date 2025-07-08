#!/usr/bin/env python3
"""Download API tools for MCP server - converts download APIs to MCP tools."""

import os
import json
import uuid
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from pydantic import BaseModel, Field, ConfigDict
from pathlib import Path
from typing import Literal, List, Dict, Optional, Any
from urllib.error import HTTPError
from eliot import start_action

DB_LITERAL = Literal[
    "pubmed", "protein", "nuccore", "ipg", "nucleotide", "structure", "genome",
    "annotinfo", "assembly", "bioproject", "biosample", "blastdbinfo", "books",
    "cdd", "clinvar", "gap", "gapplus", "grasp", "dbvar", "gene", "gds",
    "geoprofiles", "medgen", "mesh", "nlmcatalog", "omim", "orgtrack", "pmc",
    "proteinclusters", "pcassay", "protfam", "pccompound", "pcsubstance",
    "seqannot", "snp", "sra", "taxonomy", "biocollections", "gtr"
]

class EntrezDownloadRequest(BaseModel):
    ids: List[str]
    db: DB_LITERAL
    reftype: Literal["fasta", "gb"]

    model_config = {
        "json_schema_extra": {
            "example": {
                "ids": ["NM_000546.6"],
                "db": "nucleotide",
                "reftype": "fasta",
            }
        }
    }

# Type hint for local file results
LocalFileResult = Dict[Literal["path", "format", "success", "error"], Any]

class DownloadTools:
    """Handler for download-related MCP tools."""
    
    def __init__(self, mcp_server, prefix: str = "", output_dir: Optional[str] = None):
        self.mcp_server = mcp_server
        self.prefix = prefix
        self.output_dir = Path(output_dir) if output_dir else Path.cwd() / "biothings_output"
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _save_to_local_file(
        self, 
        data: Any, 
        format_type: str, 
        output_path: Optional[str] = None,
        default_prefix: str = "biothings_output"
    ) -> LocalFileResult:
        """Helper function to save data to local files.
        
        Args:
            data: The data to save
            format_type: File format ('fasta', 'gb', 'json', 'txt', etc.)
            output_path: Full output path (absolute or relative) or None to auto-generate
            default_prefix: Prefix for auto-generated filenames
            
        Returns:
            LocalFileResult: Contains path, format, success status, and optional error information
        """
        # Map format types to file extensions
        format_extensions = {
            'fasta': '.fasta',
            'gb': '.gb',
            'json': '.json',
            'txt': '.txt',
            'tsv': '.tsv',
            'alignment': '.aln'
        }
        
        extension = format_extensions.get(format_type, '.txt')
        
        if output_path is None:
            # Generate a unique filename in the default output directory
            base_name = f"{default_prefix}_{str(uuid.uuid4())[:8]}"
            file_path = self.output_dir / f"{base_name}{extension}"
        else:
            # Use the provided path
            path_obj = Path(output_path)
            if path_obj.is_absolute():
                # Absolute path - use as is, but ensure it has the right extension
                if path_obj.suffix != extension:
                    file_path = path_obj.with_suffix(extension)
                else:
                    file_path = path_obj
            else:
                # Relative path - concatenate with output directory
                if not str(output_path).endswith(extension):
                    file_path = self.output_dir / f"{output_path}{extension}"
                else:
                    file_path = self.output_dir / output_path
        
        try:
            if format_type in ['fasta', 'gb']:
                # Write sequence data
                with open(file_path, 'w') as f:
                    f.write(str(data))
            elif format_type == 'json':
                with open(file_path, 'w') as f:
                    json.dump(data, f, indent=2, default=str)
            elif format_type == 'alignment':
                # Write alignment data
                with open(file_path, 'w') as f:
                    f.write(str(data))
            else:
                # Default to text format
                with open(file_path, 'w') as f:
                    if isinstance(data, dict):
                        json.dump(data, f, indent=2, default=str)
                    else:
                        f.write(str(data))
                        
            return {
                "path": str(file_path),
                "format": format_type,
                "success": True
            }
        except Exception as e:
            return {
                "path": None,
                "format": format_type,
                "success": False,
                "error": str(e)
            }
    
    def register_tools(self):
        """Register download-related MCP tools."""
        
        @self.mcp_server.tool(
            name=f"{self.prefix}download_entrez_data",
            description="""Download data from NCBI Entrez databases using Bio.Entrez.
            
            Downloads data records from specified NCBI Entrez databases. This tool is designed to be called 
            by automated agents (like LLMs) or other services.

            **Critical Configuration:**
            The server hosting this API *must* have the `ENTREZ_EMAIL` environment variable set
            to a valid email address. NCBI requires this for Entrez queries to monitor usage
            and prevent abuse. Without it, NCBI may block requests.

            **Parameters:**
            - `ids` (List[str], required): A list of unique identifiers for the records to fetch
              from the specified Entrez database. Example: `["NM_000546.6", "AY123456.1"]`
            - `db` (DB_LITERAL, required): The target NCBI Entrez database.
              Common choices for sequences: 'nucleotide', 'protein'.
              Other examples: 'gene', 'pubmed', 'taxonomy'.
              Ensure the `ids` provided are appropriate for the selected `db`.
            - `reftype` (Literal["fasta", "gb"], required): The desired format for the
              downloaded data.
                - "fasta": Returns data in FASTA format.
                - "gb": Returns data in GenBank format.
              Ensure the chosen `reftype` is compatible with the selected `db`.

            **Returns:**
            On success: Returns the downloaded data as a single raw string with the
            data fetched from Entrez in the specified `reftype`.
            
            **Example Usage:**
            To fetch the FASTA sequence for human TP53 mRNA (NM_000546.6):
            ```
            download_entrez_data(
                ids=["NM_000546.6"],
                db="nucleotide",
                reftype="fasta"
            )
            ```
            """)
        def download_entrez_data(ids: List[str], db: DB_LITERAL, reftype: Literal["fasta", "gb"]) -> str:
            with start_action(action_type="download_entrez_data", ids=ids, db=db, reftype=reftype) as action:
                try:
                    downloaded_content = get_entrez(ids=ids, db=db, reftype=reftype)
                    action.add_success_fields(content_length=len(downloaded_content))
                    return downloaded_content
                except HTTPError as he:
                    action.add_error_fields(error=f"NCBI Entrez Error ({he.code}): {he.reason}")
                    raise ValueError(f"NCBI Entrez Error ({he.code}): {he.reason}") from he
                except Exception as e:
                    action.add_error_fields(error=f"Unexpected error during Entrez download: {str(e)}")
                    raise ValueError(f"An unexpected error occurred during Entrez download: {str(e)}") from e
        
        @self.mcp_server.tool(
            name=f"{self.prefix}download_entrez_data_local",
            description="""Download data from NCBI Entrez databases and save to local file.
            
            Same as download_entrez_data but saves the result to a local file instead of returning the content.
            This is useful for large downloads or when you want to persist the data.

            **Parameters:**
            - `ids` (List[str], required): A list of unique identifiers for the records to fetch
            - `db` (DB_LITERAL, required): The target NCBI Entrez database
            - `reftype` (Literal["fasta", "gb"], required): The desired format for the downloaded data
            - `output_path` (Optional[str]): Custom output path. If None, generates unique filename
            
            **Returns:**
            LocalFileResult containing:
            - `path`: Path to the saved file
            - `format`: File format used
            - `success`: Whether the operation succeeded
            - `error`: Error message if failed
            
            **Example Usage:**
            ```
            download_entrez_data_local(
                ids=["NM_000546.6"],
                db="nucleotide",
                reftype="fasta",
                output_path="tp53_sequence.fasta"
            )
            ```
            """)
        def download_entrez_data_local(
            ids: List[str], 
            db: DB_LITERAL, 
            reftype: Literal["fasta", "gb"],
            output_path: Optional[str] = None
        ) -> LocalFileResult:
            with start_action(action_type="download_entrez_data_local", ids=ids, db=db, reftype=reftype) as action:
                try:
                    downloaded_content = get_entrez(ids=ids, db=db, reftype=reftype)
                    result = self._save_to_local_file(
                        data=downloaded_content,
                        format_type=reftype,
                        output_path=output_path,
                        default_prefix=f"entrez_{db}"
                    )
                    action.add_success_fields(
                        content_length=len(downloaded_content),
                        saved_to=result.get("path"),
                        success=result.get("success")
                    )
                    return result
                except HTTPError as he:
                    action.add_error_fields(error=f"NCBI Entrez Error ({he.code}): {he.reason}")
                    return {
                        "path": None,
                        "format": reftype,
                        "success": False,
                        "error": f"NCBI Entrez Error ({he.code}): {he.reason}"
                    }
                except Exception as e:
                    action.add_error_fields(error=f"Unexpected error during Entrez download: {str(e)}")
                    return {
                        "path": None,
                        "format": reftype,
                        "success": False,
                        "error": f"An unexpected error occurred during Entrez download: {str(e)}"
                    }
      

def get_entrez(ids: List[str], db: DB_LITERAL, reftype: Literal["fasta", "gb"]) -> str:
    """
    Downloads data from NCBI Entrez databases.

    This function uses Bio.Entrez to fetch data based on a list of IDs from a specified database.
    
    Args:
        ids: List of unique identifiers for the records to fetch
        db: The target NCBI Entrez database
        reftype: The desired format for the downloaded data ("fasta" or "gb")
        
    Returns:
        str: The downloaded data as a string
        
    Raises:
        HTTPError: If NCBI returns an error
        Exception: For other unexpected errors
    """
    # Ensure ENTREZ_EMAIL is set
    email = os.getenv("ENTREZ_EMAIL")
    if not email:
        raise ValueError("ENTREZ_EMAIL environment variable must be set for NCBI Entrez queries")
    
    Entrez.email = email
    
    try:
        # Fetch the data
        handle = Entrez.efetch(db=db, id=ids, rettype=reftype, retmode="text")
        data = handle.read()
        handle.close()
        return data
    except HTTPError as he:
        # Re-raise HTTPError to be caught by the calling function
        raise he
    except Exception as e:
        # Re-raise other exceptions
        raise e

