"""Google Search via SerpAPI CLI"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path to import from src
sys.path.insert(0, str(Path(__file__).parent.parent))

from src import SerpAPISearch


def main():
    parser = argparse.ArgumentParser(description="Google Search using SerpAPI")
    parser.add_argument('-q', '--query', type=str, help='Search query string')
    parser.add_argument('-n', '--num_results', type=int, default=5, help='Number of search results to return (between 1 and 100)')

    args = parser.parse_args()

    try:
        search_engine = SerpAPISearch()
        results = search_engine.search(args.query, args.num_results)

        if results:
            for i, result in enumerate(results, 1):
                print(f"Result {i}")
                print(f"Title: {result.get('title', 'No title')}")
                print(f"Link: {result.get('link', 'No link')}")
                print(f"Summary: {result.get('snippet', 'No description available')}\n")
                print("---")
        else:
            print("No results found. Try a different search query.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
