#!/usr/bin/env python3
"""Download API tools for MCP server - converts download APIs to MCP tools."""

import os
import json
import uuid
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
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

class PairwiseAlignmentRequest(BaseModel):
    sequence1: str = Field(..., description="First sequence for alignment.")
    sequence2: str = Field(..., description="Second sequence for alignment.")
    match_score: float = Field(1.0, description="Score for a match.")
    mismatch_penalty: float = Field(-1.0, description="Penalty for a mismatch. Should be negative or zero.")
    open_gap_penalty: float = Field(-0.5, description="Penalty for opening a gap. Should be negative or zero.")
    extend_gap_penalty: float = Field(-0.1, description="Penalty for extending a gap. Should be negative or zero.")
    mode: Literal["global", "local"] = Field("global", description="Alignment mode: 'global' or 'local'.")

    model_config = ConfigDict(
        json_schema_extra = {
            "example": {
                "sequence1": "GATTACA",
                "sequence2": "GCATGCU",
                "match_score": 2.0,
                "mismatch_penalty": -1.0,
                "open_gap_penalty": -1.0,
                "extend_gap_penalty": -0.5,
                "mode": "global"
            }
        }
    )

class PairwiseAlignmentResponse(BaseModel):
    score: float
    aligned_sequence1: str
    aligned_sequence2: str
    full_alignment_str: str
    parameters_used: Dict

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
        
        @self.mcp_server.tool(
            name=f"{self.prefix}perform_pairwise_alignment",
            description="""Perform pairwise sequence alignment using Biopython's PairwiseAligner.
            
            Performs a pairwise sequence alignment (global or local) using Biopython's PairwiseAligner.
            You can specify sequences and alignment scoring parameters.
            
            **Parameters:**
            - `sequence1` (str): First sequence for alignment
            - `sequence2` (str): Second sequence for alignment
            - `match_score` (float): Score for a match (default: 1.0)
            - `mismatch_penalty` (float): Penalty for a mismatch (default: -1.0, should be negative or zero)
            - `open_gap_penalty` (float): Penalty for opening a gap (default: -0.5, should be negative or zero)
            - `extend_gap_penalty` (float): Penalty for extending a gap (default: -0.1, should be negative or zero)
            - `mode` (str): Alignment mode - 'global' or 'local' (default: 'global')
            
            **Returns:**
            PairwiseAlignmentResponse containing:
            - `score`: Alignment score
            - `aligned_sequence1`: First sequence with gaps
            - `aligned_sequence2`: Second sequence with gaps
            - `full_alignment_str`: Complete alignment visualization
            - `parameters_used`: Parameters used for the alignment
            
            **Example Usage:**
            ```
            perform_pairwise_alignment(
                sequence1="GATTACA",
                sequence2="GCATGCU",
                match_score=2.0,
                mismatch_penalty=-1.0,
                mode="global"
            )
            ```
            """)
        def perform_pairwise_alignment(
            sequence1: str,
            sequence2: str,
            match_score: float = 1.0,
            mismatch_penalty: float = -1.0,
            open_gap_penalty: float = -0.5,
            extend_gap_penalty: float = -0.1,
            mode: Literal["global", "local"] = "global"
        ) -> PairwiseAlignmentResponse:
            with start_action(action_type="perform_pairwise_alignment", 
                            sequence1_length=len(sequence1), 
                            sequence2_length=len(sequence2), 
                            mode=mode) as action:
                try:
                    request = PairwiseAlignmentRequest(
                        sequence1=sequence1,
                        sequence2=sequence2,
                        match_score=match_score,
                        mismatch_penalty=mismatch_penalty,
                        open_gap_penalty=open_gap_penalty,
                        extend_gap_penalty=extend_gap_penalty,
                        mode=mode
                    )
                    response = run_pairwise_alignment(request)
                    action.add_success_fields(alignment_score=response.score)
                    return response
                except ValueError as ve:
                    action.add_error_fields(error=str(ve))
                    raise ValueError(str(ve)) from ve
                except Exception as e:
                    action.add_error_fields(error=f"Unexpected error during alignment: {str(e)}")
                    raise ValueError(f"An unexpected error occurred during alignment: {str(e)}") from e
        
        @self.mcp_server.tool(
            name=f"{self.prefix}perform_pairwise_alignment_local",
            description="""Perform pairwise sequence alignment and save results to local file.
            
            Same as perform_pairwise_alignment but saves the alignment result to a local file instead of returning it.
            This is useful for preserving alignment results for further analysis.
            
            **Parameters:**
            - `sequence1` (str): First sequence for alignment
            - `sequence2` (str): Second sequence for alignment
            - `match_score` (float): Score for a match (default: 1.0)
            - `mismatch_penalty` (float): Penalty for a mismatch (default: -1.0, should be negative or zero)
            - `open_gap_penalty` (float): Penalty for opening a gap (default: -0.5, should be negative or zero)
            - `extend_gap_penalty` (float): Penalty for extending a gap (default: -0.1, should be negative or zero)
            - `mode` (str): Alignment mode - 'global' or 'local' (default: 'global')
            - `output_path` (Optional[str]): Custom output path. If None, generates unique filename
            
            **Returns:**
            LocalFileResult containing:
            - `path`: Path to the saved alignment file
            - `format`: File format used ('alignment')
            - `success`: Whether the operation succeeded
            - `error`: Error message if failed
            
            **Example Usage:**
            ```
            perform_pairwise_alignment_local(
                sequence1="GATTACA",
                sequence2="GCATGCU",
                match_score=2.0,
                mismatch_penalty=-1.0,
                mode="global",
                output_path="alignment_result.aln"
            )
            ```
            """)
        def perform_pairwise_alignment_local(
            sequence1: str,
            sequence2: str,
            match_score: float = 1.0,
            mismatch_penalty: float = -1.0,
            open_gap_penalty: float = -0.5,
            extend_gap_penalty: float = -0.1,
            mode: Literal["global", "local"] = "global",
            output_path: Optional[str] = None
        ) -> LocalFileResult:
            with start_action(action_type="perform_pairwise_alignment_local", 
                            sequence1_length=len(sequence1), 
                            sequence2_length=len(sequence2), 
                            mode=mode) as action:
                try:
                    request = PairwiseAlignmentRequest(
                        sequence1=sequence1,
                        sequence2=sequence2,
                        match_score=match_score,
                        mismatch_penalty=mismatch_penalty,
                        open_gap_penalty=open_gap_penalty,
                        extend_gap_penalty=extend_gap_penalty,
                        mode=mode
                    )
                    response = run_pairwise_alignment(request)
                    
                    # Create detailed alignment output
                    alignment_content = f"""Pairwise Alignment Results
=============================

Parameters:
- Match Score: {response.parameters_used['match_score']}
- Mismatch Penalty: {response.parameters_used['mismatch_penalty']}
- Open Gap Penalty: {response.parameters_used['open_gap_penalty']}
- Extend Gap Penalty: {response.parameters_used['extend_gap_penalty']}
- Mode: {response.parameters_used['mode']}

Alignment Score: {response.score}

Aligned Sequences:
{response.full_alignment_str}

Sequence 1 (aligned): {response.aligned_sequence1}
Sequence 2 (aligned): {response.aligned_sequence2}
"""
                    
                    result = self._save_to_local_file(
                        data=alignment_content,
                        format_type="alignment",
                        output_path=output_path,
                        default_prefix="pairwise_alignment"
                    )
                    action.add_success_fields(
                        alignment_score=response.score,
                        saved_to=result.get("path"),
                        success=result.get("success")
                    )
                    return result
                except ValueError as ve:
                    action.add_error_fields(error=str(ve))
                    return {
                        "path": None,
                        "format": "alignment",
                        "success": False,
                        "error": str(ve)
                    }
                except Exception as e:
                    action.add_error_fields(error=f"Unexpected error during alignment: {str(e)}")
                    return {
                        "path": None,
                        "format": "alignment",
                        "success": False,
                        "error": f"An unexpected error occurred during alignment: {str(e)}"
                    }

def run_pairwise_alignment(request: PairwiseAlignmentRequest) -> PairwiseAlignmentResponse:
    """Run pairwise alignment using Biopython."""
    aligner = PairwiseAligner()
    aligner.match_score = request.match_score
    aligner.mismatch_score = request.mismatch_penalty
    aligner.open_gap_score = request.open_gap_penalty
    aligner.extend_gap_score = request.extend_gap_penalty
    aligner.mode = request.mode

    seq1 = Seq(request.sequence1)
    seq2 = Seq(request.sequence2)

    alignments = list(aligner.align(seq1, seq2))

    if not alignments:
        raise ValueError("No alignment could be produced with the given sequences and parameters.")

    best_alignment = alignments[0]
    
    aligned_seq1_str = str(best_alignment[0])
    aligned_seq2_str = str(best_alignment[1])

    return PairwiseAlignmentResponse(
        score=best_alignment.score,
        aligned_sequence1=aligned_seq1_str,
        aligned_sequence2=aligned_seq2_str,
        full_alignment_str=str(best_alignment),
        parameters_used=request.model_dump()
    )

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

