import streamlit as st
from search_engines import TavilySearch


@st.cache_resource
def get_search_engine():
    return TavilySearch()


st.set_page_config(page_title="Tavily AI-Powered Search", page_icon="🔍", layout="wide")

st.sidebar.title("Tavily")
num_results = st.sidebar.slider("Number of results:", min_value=1, max_value=100, value=10)

search_query = st.text_input("Enter your search query:")

if search_query:
    st.write(f"Searching for: {search_query}")

    try:
        search_engine = get_search_engine()
        results = search_engine.search(search_query, num_results)

        if results:
            for i, result in enumerate(results, 1):
                st.subheader(f"Result {i}")
                st.write(f"**Link:** {result.get('url', 'N/A')}")
                st.write(f"**Summary:** {result.get('content', 'N/A')}")
                st.divider()
        else:
            st.write("No results found.")
    except Exception as e:
        st.error(f"Error: {e}")
