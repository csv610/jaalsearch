import requests
import os
import argparse

# Set your SerpAPI key as an environment variable or directly here
SERPAPI_API_KEY = os.getenv('SERPAPI_API_KEY')

def search_google(query, num_results=5):
    """
    Perform a Google search using SerpAPI.
    
    Args:
        query (str): Search query string.
        num_results (int): Number of search results to return (between 1 and 100).
        
    Returns:
        list: A list of search results dictionaries.
    """
    # Input validation
    if not query.strip():
        print("Error: Search query cannot be empty.")
        return []
    if not (1 <= num_results <= 100):
        print("Error: Number of results must be between 1 and 100.")
        return []
    
    try:
        # Fetching results from SerpAPI using Python API
        params = {
            "q": query,
            "num": num_results,
            "api_key": SERPAPI_API_KEY,
            "engine": "google"
        }
        response = requests.get('https://serpapi.com/search', params=params).json()
        results = response.get('organic_results', [])
        return results
    except Exception as e:
        # Handle any potential errors in the search process
        print(f"An error occurred: {e}")
        return []

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Google Search using SerpAPI")
    parser.add_argument('-q', '--query', type=str, help='Search query string')
    parser.add_argument('-n', '--num_results', type=int, default=5, help='Number of search results to return (between 1 and 100)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Get search results
    results = search_google(args.query, args.num_results)
    
    # Display results if any are found
    if results:
        for i, result in enumerate(results, 1):
            print(f"Result {i}")
            print(f"Title: {result.get('title', 'No title')}")
            print(f"Link: {result.get('link', 'No link')}")
            print(f"Summary: {result.get('snippet', 'No description available')}\n")
            print("---")
    else:
        print("No results found. Try a different search query.")

if __name__ == "__main__":
    main()
