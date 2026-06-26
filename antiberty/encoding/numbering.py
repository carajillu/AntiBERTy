
import re
from typing import Optional, Dict, Union


def extract_variable_domain_anarci(
    sequence: str,
    scheme: str = "imgt",
    chain_type: Optional[str] = None,
    return_details: bool = False,
) -> Union[str, Dict[str, object]]:
    """
    Extract an antibody variable domain from a full heavy or light chain sequence
    using ANARCI.

    Parameters
    ----------
    sequence : str
        Full antibody amino-acid sequence. May include signal peptide and/or
        constant region.
    scheme : str
        ANARCI numbering scheme. Common options include:
        "imgt", "chothia", "kabat", "martin", "aho", "wolfguy".
    chain_type : str, optional
        Restrict chain type recognition.
        Options:
            "heavy" or "H"
            "kappa" or "K"
            "lambda" or "L"
            "light"  -> allows both K and L
        If None, allows H, K, and L.
    return_details : bool
        If True, return a dictionary with sequence, chain type, scheme,
        numbering, and residue count.

    Returns
    -------
    str or dict
        Variable-domain sequence, or details dictionary if return_details=True.

    Requirements
    ------------
    pip install ANARCI

    or, depending on your environment, install ANARCI with its HMMER dependency
    via conda/bioconda.
    """

    try:
        from anarci import number
    except ImportError as exc:
        raise ImportError(
            "ANARCI is not installed. Install it with `pip install ANARCI` "
            "or via conda/bioconda, depending on your environment."
        ) from exc

    if not sequence:
        raise ValueError("Input sequence is empty.")

    # Remove FASTA header, whitespace, numbering, and non-amino-acid characters.
    lines = sequence.strip().splitlines()
    lines = [line for line in lines if not line.startswith(">")]
    clean_seq = "".join(lines).upper()
    clean_seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYXBZUO]", "", clean_seq)

    if len(clean_seq) < 70:
        raise ValueError(
            "Sequence is probably too short to contain a complete antibody variable domain."
        )

    # Map user-friendly chain type names to ANARCI chain codes.
    # H = heavy, K = kappa, L = lambda.
    if chain_type is None:
        allow = {"H", "K", "L"}
    else:
        ct = chain_type.lower()
        if ct in {"heavy", "h"}:
            allow = {"H"}
        elif ct in {"kappa", "k"}:
            allow = {"K"}
        elif ct in {"lambda", "l"}:
            allow = {"L"}
        elif ct == "light":
            allow = {"K", "L"}
        else:
            raise ValueError(
                "chain_type must be one of None, 'heavy', 'kappa', "
                "'lambda', 'light', 'H', 'K', or 'L'."
            )

    numbering, inferred_chain_type = number(
        clean_seq,
        scheme=scheme,
        allow=allow,
    )

    if not numbering:
        raise ValueError(
            "ANARCI could not identify a valid antibody variable domain "
            "in this sequence."
        )

    # ANARCI numbering is a list of entries like:
    # [((position_number, insertion_code), amino_acid), ...]
    # Gaps are normally represented by '-'. We remove them.
    variable_domain = "".join(
        aa for position, aa in numbering
        if aa not in {"-", "."}
    )

    if not variable_domain:
        raise ValueError(
            "ANARCI returned numbering, but no amino-acid residues could be extracted."
        )

    if return_details:
        return {
            "variable_domain": variable_domain,
            "length": len(variable_domain),
            "chain_type": inferred_chain_type,
            "scheme": scheme,
            "numbering": numbering,
            "input_length": len(clean_seq),
        }

    return variable_domain

import re
from typing import Optional, Dict, Union


def pad_variable_domain_anarci(
    variable_domain: str,
    target_length: int,
    scheme: str = "imgt",
    chain_type: Optional[str] = None,
    pad_char: str = "-",
    pad_position: str = "right",
    if_longer: str = "raise",
    return_details: bool = False,
) -> Union[str, Dict[str, object]]:
    """
    Number an antibody variable domain with ANARCI and pad it to a fixed length.

    Parameters
    ----------
    variable_domain : str
        Antibody variable-domain amino-acid sequence.
        This should already be the VH, VK, or VL sequence rather than the full chain.

    target_length : int
        Desired final sequence length after ANARCI alignment and additional padding.

    scheme : str
        ANARCI numbering scheme.
        Common options include:
            "imgt", "chothia", "kabat", "martin", "aho", "wolfguy"

    chain_type : str, optional
        Restrict chain type recognition.
        Options:
            None      -> allow H, K, and L
            "heavy"   or "H"
            "kappa"   or "K"
            "lambda"  or "L"
            "light"   -> allow K and L

    pad_char : str
        Character used for padding. Default is "-".

    pad_position : str
        Where to add extra padding if the ANARCI-aligned sequence is shorter
        than target_length.
        Options:
            "right"  -> append padding to the C-terminus
            "left"   -> prepend padding to the N-terminus
            "both"   -> split padding between both ends

    if_longer : str
        What to do if the ANARCI-aligned sequence is longer than target_length.
        Options:
            "raise"    -> raise ValueError
            "truncate" -> truncate to target_length

    return_details : bool
        If True, return a dictionary containing the padded sequence and metadata.

    Returns
    -------
    str or dict
        Padded variable-domain sequence, or details dictionary if return_details=True.

    Requirements
    ------------
    ANARCI must be installed and importable:

        from anarci import number

    Notes
    -----
    This function preserves ANARCI-introduced internal gaps, then adds extra
    padding only if needed to reach target_length.
    """

    try:
        from anarci import number
    except ImportError as exc:
        raise ImportError(
            "ANARCI is not installed or not importable. Install ANARCI and "
            "its dependencies, including HMMER, before using this function."
        ) from exc

    if not variable_domain:
        raise ValueError("Input variable_domain sequence is empty.")

    if target_length <= 0:
        raise ValueError("target_length must be a positive integer.")

    if len(pad_char) != 1:
        raise ValueError("pad_char must be a single character.")

    pad_position = pad_position.lower()
    if pad_position not in {"left", "right", "both"}:
        raise ValueError("pad_position must be one of: 'left', 'right', or 'both'.")

    if_longer = if_longer.lower()
    if if_longer not in {"raise", "truncate"}:
        raise ValueError("if_longer must be one of: 'raise' or 'truncate'.")

    # Clean input sequence.
    seq = variable_domain.strip().upper()
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYXBZUO]", "", seq)

    if len(seq) < 70:
        raise ValueError(
            "Input sequence looks too short for a complete antibody variable domain."
        )

    # Map user-friendly chain type to ANARCI chain codes.
    if chain_type is None:
        allow = {"H", "K", "L"}
    else:
        ct = chain_type.lower()
        if ct in {"heavy", "h"}:
            allow = {"H"}
        elif ct in {"kappa", "k"}:
            allow = {"K"}
        elif ct in {"lambda", "l"}:
            allow = {"L"}
        elif ct == "light":
            allow = {"K", "L"}
        else:
            raise ValueError(
                "chain_type must be one of None, 'heavy', 'kappa', "
                "'lambda', 'light', 'H', 'K', or 'L'."
            )

    # Run ANARCI numbering.
    numbering, inferred_chain_type = number(
        seq,
        scheme=scheme,
        allow=allow,
    )

    if not numbering:
        raise ValueError(
            "ANARCI could not number this sequence as an antibody variable domain."
        )

    # Reconstruct ANARCI-aligned sequence.
    # ANARCI gap positions are usually represented as '-'.
    anarci_aligned = "".join(
        aa if aa not in {".", " "} else pad_char
        for position, aa in numbering
    )

    # Normalize any ANARCI gaps to the requested pad character.
    anarci_aligned = anarci_aligned.replace("-", pad_char)

    aligned_length = len(anarci_aligned)

    # Handle sequences longer than target_length.
    if aligned_length > target_length:
        if if_longer == "raise":
            raise ValueError(
                f"ANARCI-aligned sequence length is {aligned_length}, "
                f"which is longer than target_length={target_length}."
            )
        elif if_longer == "truncate":
            padded = anarci_aligned[:target_length]

    # Handle sequences equal to target_length.
    elif aligned_length == target_length:
        padded = anarci_aligned

    # Add extra padding if needed.
    else:
        n_pad = target_length - aligned_length

        if pad_position == "right":
            padded = anarci_aligned + pad_char * n_pad

        elif pad_position == "left":
            padded = pad_char * n_pad + anarci_aligned

        elif pad_position == "both":
            left_pad = n_pad // 2
            right_pad = n_pad - left_pad
            padded = pad_char * left_pad + anarci_aligned + pad_char * right_pad

    if return_details:
        return {
            "padded_sequence": padded,
            "original_sequence": seq,
            "original_length": len(seq),
            "anarci_aligned_sequence": anarci_aligned,
            "anarci_aligned_length": aligned_length,
            "target_length": target_length,
            "scheme": scheme,
            "chain_type": inferred_chain_type,
            "numbering": numbering,
        }

    return padded