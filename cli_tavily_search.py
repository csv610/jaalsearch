import argparse
from search_engines import TavilySearch


def main():
    parser = argparse.ArgumentParser(description="Tavily AI-Powered Search CLI")
    parser.add_argument('-q', '--query', type=str, help='The search query string')
    parser.add_argument('-n', '--num-results', type=int, default=5, help='Number of results to fetch')
    args = parser.parse_args()

    try:
        search_engine = TavilySearch()
        results = search_engine.search(args.query, args.num_results)

        if results:
            for i, result in enumerate(results, 1):
                print(f"Result {i}")
                print(f"Link: {result.get('url', 'N/A')}")
                print(f"Summary: {result.get('content', 'N/A')}\n")
                print("-" * 40)
        else:
            print("No results found.")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
