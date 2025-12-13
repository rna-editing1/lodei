import pysam

def is_bam_file(filepath):
    """
    Check if the given filepath is a valid BAM file.
    
    Args:
        filepath (str): Path to the file to check
        
    Returns:
        bool: True if the file is a valid BAM file, False otherwise
    """
    try:
        con = pysam.AlignmentFile(filepath, "rb")
        con.close()
        return True
    except Exception:
        return False

def is_supported_librarytype(libtype):
    """
    Check if the given library type is supported.
    
    Args:
        libtype (str): Library type to check
    Returns:
        bool: True if supported, False otherwise
    """
    supported_types = ["U", "SF", "SR", "ISR", "ISF"]
    return libtype in supported_types