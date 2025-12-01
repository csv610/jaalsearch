"""DuckDuckGo Search CLI"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path to import from src
sys.path.insert(0, str(Path(__file__).parent.parent))

from src import DuckDuckGoSearch


def main():
    parser = argparse.ArgumentParser(description="DuckDuckGo Search CLI")
    parser.add_argument("-q", "--query", type=str, help="The search query")
    parser.add_argument("-n", "--num_results", type=int, default=5, help="Number of results to display (default: 5)")
    parser.add_argument("-s", "--safesearch", type=str, choices=["Off", "Moderate", "Strict"], default="Moderate", help="SafeSearch level (default: Moderate)")
    args = parser.parse_args()

    try:
        search_engine = DuckDuckGoSearch()
        results = search_engine.search(args.query, args.num_results, safesearch=args.safesearch)

        if results:
            for i, result in enumerate(results, 1):
                print(f"Result {i}: {result['title']}")
                print(f"Link: {result['href']}")
                print(f"Summary: {result['body']}\n")
        else:
            print("No results found.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
