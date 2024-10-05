import argparse
import os
from tavily import TavilyClient

# Ensure you have set the TAVILY_API_KEY environment variable
tavily_api_key = os.getenv("TAVILY_API_KEY")
if not tavily_api_key:
    raise EnvironmentError("Please set the TAVILY_API_KEY environment variable.")

def get_tavily_client():
    return TavilyClient(tavily_api_key)

def search_tavily(client, query, num_results=5):
    try:
        response = client.search(query, max_results=num_results)
        return response.get('results', [])
    except Exception as e:
        print(f"An error occurred while fetching search results: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Tavily AI-Powered Search CLI")
    parser.add_argument('-q', '--query', type=str, help='The search query string')
    parser.add_argument('-n', '--num-results', type=int, default=5, help='Number of results to fetch')
    args = parser.parse_args()

    tavily_client = get_tavily_client()
    print(f"Searching for: {args.query}\n")
    results = search_tavily(tavily_client, args.query, args.num_results)

    if results:
        for i, result in enumerate(results, 1):
            print(f"Result {i}")
            print(f"Link: {result.get('url', 'N/A')}")
            print(f"Summary: {result.get('content', 'N/A')}\n")
            print("-" * 40)
    else:
        print("No results found.")

if __name__ == "__main__":
    main()
