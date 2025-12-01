"""
Helper functions for JAAlSearch.

Utility functions for common operations across the application.
"""

from typing import List, Dict, Any, Optional


def format_results(results: List[Dict[str, Any]], format_type: str = "text") -> str:
    """
    Format search results for display.

    Args:
        results: List of result dictionaries.
        format_type: Format type ('text', 'json', 'csv').

    Returns:
        Formatted results string.
    """
    if format_type == "json":
        import json
        return json.dumps(results, indent=2)
    elif format_type == "csv":
        return format_as_csv(results)
    else:
        return format_as_text(results)


def format_as_text(results: List[Dict[str, Any]]) -> str:
    """Format results as plain text."""
    output = []
    for i, result in enumerate(results, 1):
        output.append(f"Result {i}:")
        for key, value in result.items():
            output.append(f"  {key}: {value}")
        output.append("")
    return "\n".join(output)


def format_as_csv(results: List[Dict[str, Any]]) -> str:
    """Format results as CSV."""
    if not results:
        return ""

    import csv
    from io import StringIO

    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)
    return output.getvalue()


def validate_query(query: str, min_length: int = 1, max_length: int = 1000) -> bool:
    """
    Validate search query.

    Args:
        query: Search query string.
        min_length: Minimum query length.
        max_length: Maximum query length.

    Returns:
        True if query is valid, False otherwise.
    """
    if not query or not isinstance(query, str):
        return False

    query = query.strip()
    if len(query) < min_length or len(query) > max_length:
        return False

    return True


def truncate_text(text: str, max_length: int = 100, suffix: str = "...") -> str:
    """
    Truncate text to maximum length.

    Args:
        text: Text to truncate.
        max_length: Maximum length.
        suffix: Suffix to add if truncated.

    Returns:
        Truncated text.
    """
    if len(text) <= max_length:
        return text

    return text[:max_length - len(suffix)] + suffix


def remove_duplicates(results: List[Dict[str, Any]], key: str = "title") -> List[Dict[str, Any]]:
    """
    Remove duplicate results based on key field.

    Args:
        results: List of result dictionaries.
        key: Key field to check for duplicates.

    Returns:
        List with duplicates removed.
    """
    seen = set()
    unique_results = []

    for result in results:
        value = result.get(key, "")
        if value not in seen:
            seen.add(value)
            unique_results.append(result)

    return unique_results


def filter_results(
    results: List[Dict[str, Any]],
    filters: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    Filter results based on criteria.

    Args:
        results: List of result dictionaries.
        filters: Dictionary of filter criteria.

    Returns:
        Filtered results.
    """
    filtered = results

    for field, pattern in filters.items():
        filtered = [
            r for r in filtered
            if pattern.lower() in str(r.get(field, "")).lower()
        ]

    return filtered


def get_result_summary(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Get summary statistics about results.

    Args:
        results: List of result dictionaries.

    Returns:
        Summary dictionary with count and other stats.
    """
    return {
        "total_results": len(results),
        "has_results": len(results) > 0,
        "fields": list(results[0].keys()) if results else [],
    }
