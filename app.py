import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import re

# PAGE CONFIGURATION
st.set_page_config(
    page_title="GeneScope - GC Analyzer",
    page_icon="🧬",
    layout="wide"
)

# CUSTOM WEBSITE DESIGN
st.markdown("""
<style>
.stApp {
    background-color: #f5faf8;
}
h1 {
    color: #155e52;
}
h2, h3 {
    color: #26766c;
}
</style>
""", unsafe_allow_html=True)


# FUNCTION TO VALIDATE DNA
def clean_sequence(raw):

    lines = []

    for line in raw.splitlines():
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        lines.append(re.sub(r"\s+", "", line).upper())

    sequence = "".join(lines)

    if not sequence:
        raise ValueError("Please enter a DNA sequence.")

    invalid = set(sequence) - set("ATGC")

    if invalid:
        raise ValueError("Invalid DNA bases! Use only A, T, G and C.")

    return sequence


# FUNCTION TO CALCULATE GC CONTENT
def gc_calculator(sequence):

    total = len(sequence)

    A = sequence.count("A")
    T = sequence.count("T")
    G = sequence.count("G")
    C = sequence.count("C")

    gc = ((G + C) / total) * 100
    at = ((A + T) / total) * 100

    return {
        "Length": total,
        "A": A,
        "T": T,
        "G": G,
        "C": C,
        "GC Content (%)": round(gc, 2),
        "AT Content (%)": round(at, 2)
    }


# WEBSITE HEADER
st.title("🧬 GeneScope")
st.subheader("DNA GC Content Calculator & Comparative Analysis")

st.write(
    "A Python-based bioinformatics tool for analysing "
    "DNA nucleotide composition and GC content."
)

st.divider()


# WEBSITE TABS
tab1, tab2, tab3 = st.tabs([
    "Single Sequence Analysis",
    "Comparative Analysis",
    "Sliding Window Analysis"
])


# TAB 1: SINGLE SEQUENCE
with tab1:

    st.header("DNA Sequence Analysis")

    dna_input = st.text_area(
        "Enter your DNA sequence:",
        value="ATGCGCATTA",
        height=120
    )

    if st.button("Calculate GC Content"):

        try:
            sequence = clean_sequence(dna_input)
            result = gc_calculator(sequence)

            st.success("DNA sequence analysed successfully!")

            col1, col2, col3 = st.columns(3)

            col1.metric("Sequence Length", result["Length"])
            col2.metric("GC Content", str(result["GC Content (%)"]) + "%")
            col3.metric("AT Content", str(result["AT Content (%)"]) + "%")

            st.subheader("Nucleotide Composition")

            st.write(result)

            # NUCLEOTIDE GRAPH
            bases = ["A", "T", "G", "C"]
            counts = [
                result["A"],
                result["T"],
                result["G"],
                result["C"]
            ]

            fig, ax = plt.subplots()

            ax.bar(bases, counts)

            ax.set_xlabel("Nucleotides")
            ax.set_ylabel("Count")
            ax.set_title("DNA Nucleotide Composition")

            st.pyplot(fig)
            plt.close(fig)

        except ValueError as error:
            st.error(str(error))


# TAB 2: COMPARATIVE ANALYSIS
with tab2:

    st.header("Comparative GC Content Analysis")

    sequences = {
        "Sequence 1": "ATGCGCATTA",
        "Sequence 2": "GCGCGCATGC",
        "Sequence 3": "ATATTAATTA"
    }

    results = []

    for name, dna in sequences.items():

        result = gc_calculator(dna)
        result["Sequence"] = name

        results.append(result)

    df = pd.DataFrame(results)

    st.subheader("Analysis Results")

    st.dataframe(df, use_container_width=True)

    # COMPARISON GRAPH
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.bar(df["Sequence"], df["GC Content (%)"])

    ax.set_xlabel("DNA Sequences")
    ax.set_ylabel("GC Content (%)")
    ax.set_title("Comparative GC Content")

    ax.set_ylim(0, 100)
    ax.grid(axis="y", alpha=0.3)

    st.pyplot(fig)
    plt.close(fig)

    highest = df.loc[df["GC Content (%)"].idxmax()]
    lowest = df.loc[df["GC Content (%)"].idxmin()]

    st.write(
        "Highest GC Content:",
        highest["Sequence"],
        highest["GC Content (%)"], "%"
    )

    st.write(
        "Lowest GC Content:",
        lowest["Sequence"],
        lowest["GC Content (%)"], "%"
    )

    st.download_button(
        "Download Results",
        df.to_csv(index=False),
        file_name="gc_analysis.csv",
        mime="text/csv"
    )


# TAB 3: SLIDING WINDOW ANALYSIS
with tab3:

    st.header("Sliding Window GC Analysis")

    dna = st.text_area(
        "Enter a longer DNA sequence:",
        value="ATGCGCATTA" * 20
    )

    window = st.number_input(
        "Window Size",
        min_value=2,
        value=10
    )

    step = st.number_input(
        "Step Size",
        min_value=1,
        value=5
    )

    if st.button("Analyse GC Variation"):

        try:
            sequence = clean_sequence(dna)

            if window > len(sequence):
                st.error("Window size exceeds sequence length.")

            else:
                positions = []
                values = []

                for start in range(
                    0,
                    len(sequence) - window + 1,
                    step
                ):

                    part = sequence[start:start + window]

                    gc = (
                        (part.count("G") + part.count("C"))
                        / len(part)
                    ) * 100

                    positions.append(start + 1)
                    values.append(round(gc, 2))

                graph_data = pd.DataFrame({
                    "Position": positions,
                    "GC Content": values
                })

                st.dataframe(graph_data)

                fig, ax = plt.subplots()

                ax.plot(positions, values, marker="o")

                ax.set_xlabel("Sequence Position")
                ax.set_ylabel("GC Content (%)")
                ax.set_title("Sliding Window GC Analysis")

                ax.set_ylim(0, 100)
                ax.grid(alpha=0.3)

                st.pyplot(fig)
                plt.close(fig)

        except ValueError as error:
            st.error(str(error))


# FOOTER
st.divider()

st.caption(
    "GeneScope | A Student Bioinformatics Project | "
    "Developed using Python and Streamlit"
)

