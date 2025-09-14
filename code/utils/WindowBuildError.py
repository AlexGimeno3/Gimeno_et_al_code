class WindowBuildError(Exception):
    """Custom exception for window building errors."""
    def __init__(self, message):
        super().__init__(message)