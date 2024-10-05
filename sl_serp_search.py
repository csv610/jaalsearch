import streamlit as st
import requests
import os

# Set your SerpAPI key as an environment variable or directly here
SERPAPI_API_KEY = os.getenv('SERPAPI_API_KEY')

# Function to perform a Google search using SerpAPI Python API
def search_google(query, num_results=5):
    # Input validation
    if not query.strip():
        st.error("Search query cannot be empty.")
        return []
    if not (1 <= num_results <= 100):
        st.error("Number of results must be between 1 and 100.")
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
        st.error(f"An error occurred: {e}")
        return []

# Set wide layout for the app with title and icon
st.set_page_config(layout="wide", page_title="Google Search", page_icon="🔍")

# Sidebar for search settings
st.sidebar.title("SerpAPI")
num_results = st.sidebar.slider("Number of results:", min_value=1, max_value=100, value=10)

# Input field for the search query
search_query = st.text_input("Enter your search query:", on_change=lambda: st.session_state.update({'button_pressed': True}))

# Automatically search when Enter is pressed
if search_query.strip() and st.session_state.get('button_pressed'):
    st.write(f"🔍 Searching for: **{search_query}**")
    results = search_google(search_query, num_results)
    
    # Display results if any are found
    if results:
        for i, result in enumerate(results, 1):
            st.subheader(f"Result {i}")
            st.write(f"**Title:** {result.get('title', 'No title')}")
            st.write(f"**Link:** {result.get('link', 'No link')}")
            st.write(f"**Summary:** {result.get('snippet', 'No description available')}")
            st.write("---")
    else:
        st.warning("No results found. Try a different search query.")
