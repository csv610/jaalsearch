import argparse
from duckduckgo_search import DDGS

# Perform a search using DuckDuckGo
def search_duckduckgo(ddgs, query, num_results=5, safesearch="Moderate"):
    try:
        return list(ddgs.text(query, max_results=num_results, safesearch=safesearch))
    except Exception as e:
        print(f"An error occurred while searching: {e}")
        return []

if __name__ == "__main__":
    # Argument parser for CLI inputs
    parser = argparse.ArgumentParser(description="DuckDuckGo Search CLI")
    parser.add_argument("-q", "--query", type=str, help="The search query")
    parser.add_argument("-n", "--num_results", type=int, default=5, help="Number of results to display (default: 5)")
    parser.add_argument("-s", "--safesearch", type=str, choices=["Off", "Moderate", "Strict"], default="Moderate", help="SafeSearch level (default: Moderate)")
    args = parser.parse_args()

    # Create a DDGS instance
    ddgs = DDGS()

    # Perform the search
    results = search_duckduckgo(ddgs, args.query, args.num_results, safesearch=args.safesearch)

    # Display the search results
    if results:
        for i, result in enumerate(results, 1):
            print(f"Result {i}: {result['title']}")
            print(f"Link: {result['href']}")
            print(f"Summary: {result['body']}\n")
    else:
        print("No results found.")
