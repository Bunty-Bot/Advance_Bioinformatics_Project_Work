I want the proper docstring and the annotations for this code: import streamlit as st
import json
import os
import pandas as pd
import numpy as np
import altair as alt
import re

# --- CONFIGURATION ---
BASE_FASTP_DIR = "/lustre/BIF/nobackup/BIF30806/pizza_pepperoni/output_directory/fastp_results"
BASE_STRINGTIE_DIR = "/lustre/BIF/nobackup/BIF30806/pizza_pepperoni/output_directory/stringtie_results"
BASE_DESEQ2_DIR = "/lustre/BIF/nobackup/BIF30806/pizza_pepperoni/output_directory/deseq2_results"
BASE_HISAT2_DIR = "/lustre/BIF/nobackup/BIF30806/pizza_pepperoni/output_directory/hisat2_results"

# --- PAGE SETUP ---
st.set_page_config(
    page_title="Pizza Pepperoni Dashboard",
    page_icon="🍕",
    layout="wide"
)


# --- HELPER: CUSTOM CHARTS ---
def make_dynamic_line_chart(df, x_col, y_col, title, y_title, color_hex):
    if x_col not in df.columns and x_col == df.index.name:
        df = df.reset_index()
    chart = alt.Chart(df).mark_line(color=color_hex, strokeWidth=2).encode(
        x=alt.X(x_col, title="Position (bp)"),
        y=alt.Y(y_col, title=y_title, scale=alt.Scale(zero=False, padding=1)),
        tooltip=[x_col, y_col]
    ).properties(title=title, height=300).interactive()
    return chart


def make_dynamic_scatter(df, x, y, color, size, tooltip_cols):
    chart = alt.Chart(df).mark_circle().encode(
        x=alt.X(x, scale=alt.Scale(zero=False, padding=1)),
        y=alt.Y(y, scale=alt.Scale(zero=False, padding=1)),
        color=color,
        size=size,
        tooltip=tooltip_cols
    ).interactive()
    return chart


# --- DATA PARSING FUNCTIONS ---
def parse_hisat2_report(file_path):
    metrics = {}
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        total_match = re.search(r"(\d+) reads; of these:", content)
        if total_match: metrics['total_reads'] = int(total_match.group(1))
        rate_match = re.search(r"([\d\.]+)% overall alignment rate", content)
        if rate_match: metrics['alignment_rate'] = float(rate_match.group(1))
        unique_match = re.search(r"(\d+) \([\d\.]+%?\) aligned concordantly exactly 1 time", content)
        if unique_match: metrics['unique_aligned'] = int(unique_match.group(1))
        multi_match = re.search(r"(\d+) \([\d\.]+%?\) aligned concordantly >1 times", content)
        if multi_match: metrics['multi_aligned'] = int(multi_match.group(1))
        unaligned_conc_match = re.search(r"(\d+) \([\d\.]+%?\) aligned concordantly 0 times", content)
        if unaligned_conc_match: metrics['unaligned_concordant'] = int(unaligned_conc_match.group(1))
    except Exception as e:
        print(f"Error parsing HISAT2: {e}")
    return metrics


# --- DATA LOADING FUNCTIONS ---
@st.cache_data
def load_cohort_summary():
    summary_data = []
    if not os.path.exists(BASE_FASTP_DIR): return pd.DataFrame()
    sample_dirs = [d for d in os.listdir(BASE_FASTP_DIR) if os.path.isdir(os.path.join(BASE_FASTP_DIR, d))]
    sample_dirs.sort()

    for sample_id in sample_dirs:
        record = {"Sample ID": sample_id}
        # 1. FastP
        fastp_path = os.path.join(BASE_FASTP_DIR, sample_id, f"{sample_id}.json")
        if os.path.exists(fastp_path):
            try:
                with open(fastp_path, 'r') as f:
                    d = json.load(f)
                    stats = d.get('summary', {}).get('before_filtering', {})
                    record["Total Reads (M)"] = stats.get('total_reads', 0) / 1_000_000
                    record["GC Content (%)"] = stats.get('gc_content', 0) * 100
                    record["Q30 Rate (%)"] = stats.get('q30_rate', 0) * 100
            except:
                pass
        # 2. StringTie
        st_path = os.path.join(BASE_STRINGTIE_DIR, sample_id, f"{sample_id}.abundance.txt")
        if os.path.exists(st_path):
            try:
                df_st = pd.read_csv(st_path, sep='\t')
                genes_detected = df_st[df_st['TPM'] > 1].shape[0]
                record["Genes Detected (>1 TPM)"] = genes_detected
                record["Max TPM"] = df_st['TPM'].max()
                record["Mean TPM"] = df_st['TPM'].mean()
            except:
                pass
        # 3. HISAT2
        hisat_path = os.path.join(BASE_HISAT2_DIR, sample_id, f"{sample_id}_report.txt")
        if os.path.exists(hisat_path):
            h_metrics = parse_hisat2_report(hisat_path)
            if 'alignment_rate' in h_metrics:
                record["Alignment Rate (%)"] = h_metrics['alignment_rate']
        summary_data.append(record)
    return pd.DataFrame(summary_data)


def load_stringtie_detail(sample_id):
    path = os.path.join(BASE_STRINGTIE_DIR, sample_id, f"{sample_id}.abundance.txt")
    if os.path.exists(path): return pd.read_csv(path, sep='\t')
    return None


def load_fastp_detail(sample_id):
    path = os.path.join(BASE_FASTP_DIR, sample_id, f"{sample_id}.json")
    if os.path.exists(path):
        with open(path, 'r') as f: return json.load(f)
    return None


def load_hisat2_detail(sample_id):
    path = os.path.join(BASE_HISAT2_DIR, sample_id, f"{sample_id}_report.txt")
    if os.path.exists(path): return parse_hisat2_report(path)
    return None


@st.cache_data
def load_deseq2_files():
    if not os.path.exists(BASE_DESEQ2_DIR): return []
    files = [f for f in os.listdir(BASE_DESEQ2_DIR) if f.endswith(".csv")]
    files.sort()
    return files


@st.cache_data
def process_deseq2_data(filename):
    path = os.path.join(BASE_DESEQ2_DIR, filename)
    try:
        df = pd.read_csv(path)
        df = df.dropna(subset=['padj', 'log2FoldChange'])
        df['neg_log10_padj'] = -np.log10(df['padj'] + 1e-300)
        conditions = [(df['padj'] < 0.05) & (df['log2FoldChange'] > 1),
                      (df['padj'] < 0.05) & (df['log2FoldChange'] < -1)]
        choices = ['Upregulated', 'Downregulated']
        df['Regulation'] = np.select(conditions, choices, default='Not Significant')
        return df
    except Exception as e:
        return None


# --- MAIN UI ---
def main():
    st.title("🍕 Pizza Pepperoni Dashboard")
    df_summary = load_cohort_summary()
    if df_summary.empty:
        st.error(f"No data found on server.")
        return

    view_mode = st.radio("View Mode:", ["📊 Cohort Overview & Comparisons", "🔍 Single Sample Details"], horizontal=True)
    st.divider()

    # ==========================================
    # VIEW 1: COHORT OVERVIEW
    # ==========================================
    if view_mode == "📊 Cohort Overview & Comparisons":
        tab_qc, tab_align, tab_exp, tab_dea = st.tabs(
            ["🔍 Quality Control", "🧬 Alignment (HISAT2)", "🧪 Expression", "🌋 Differential Expression"])

        with tab_qc:
            st.subheader("Sequencing Quality Landscape")
            c1, c2 = st.columns(2)
            with c1:
                st.markdown("##### Sequencing Depth")
                st.bar_chart(df_summary.set_index("Sample ID")["Total Reads (M)"], color="#2563eb")
            with c2:
                st.markdown("##### Quality vs GC Content")
                if "GC Content (%)" in df_summary.columns:
                    chart = make_dynamic_scatter(
                        df_summary, x='GC Content (%)', y='Q30 Rate (%)',
                        color='Sample ID', size='Total Reads (M)',
                        tooltip_cols=['Sample ID', 'Total Reads (M)', 'GC Content (%)', 'Q30 Rate (%)']
                    )
                    st.altair_chart(chart, use_container_width=True)

        with tab_align:
            st.subheader("Alignment Performance")
            if "Alignment Rate (%)" in df_summary.columns:
                st.markdown("##### Overall Alignment Rate (%)")
                st.caption("Axis is zoomed to the data range to show differences.")

                # Dynamic Bar Chart for Alignment
                base_chart = alt.Chart(df_summary).mark_bar().encode(
                    x=alt.X("Sample ID", sort=None),
                    # UPDATED: zero=False to zoom in on differences
                    y=alt.Y("Alignment Rate (%)", scale=alt.Scale(zero=False)),
                    color=alt.condition(
                        alt.datum['Alignment Rate (%)'] < 80,
                        alt.value('red'),
                        alt.value('#10b981')
                    ),
                    tooltip=["Sample ID", "Alignment Rate (%)"]
                ).interactive()
                st.altair_chart(base_chart, use_container_width=True)
            else:
                st.warning("No HISAT2 data found.")

        with tab_exp:
            st.subheader("Transcriptome Complexity")
            if "Genes Detected (>1 TPM)" in df_summary.columns:
                c1, c2 = st.columns(2)
                with c1:
                    st.markdown("##### Detected Genes (TPM > 1)")
                    st.caption("Axis is zoomed to the data range.")
                    # UPDATED: Replaced simple st.bar_chart with custom Altair for dynamic scaling
                    chart_genes = alt.Chart(df_summary).mark_bar(color="#059669").encode(
                        x=alt.X("Sample ID", sort=None),
                        y=alt.Y("Genes Detected (>1 TPM)", scale=alt.Scale(zero=False)),
                        tooltip=["Sample ID", "Genes Detected (>1 TPM)"]
                    ).interactive()
                    st.altair_chart(chart_genes, use_container_width=True)

                with c2:
                    st.markdown("##### Expression Spread")
                    chart = make_dynamic_scatter(
                        df_summary, x="Mean TPM", y="Max TPM",
                        color="Sample ID", size="Genes Detected (>1 TPM)",
                        tooltip_cols=['Sample ID', 'Mean TPM', 'Max TPM']
                    )
                    st.altair_chart(chart, use_container_width=True)
            else:
                st.warning("No StringTie data found.")

        with tab_dea:
            st.header("Differential Expression Analysis")
            deseq_files = load_deseq2_files()
            if not deseq_files:
                st.warning(f"No CSV files found.")
            else:
                selected_file = st.selectbox("Select Comparison:", deseq_files)
                df_volcano = process_deseq2_data(selected_file)
                if df_volcano is not None:
                    c1, c2 = st.columns([3, 1])
                    with c1:
                        st.subheader("Volcano Plot")
                        volcano_chart = alt.Chart(df_volcano).mark_circle(size=60).encode(
                            x=alt.X('log2FoldChange', title='Log2 Fold Change', scale=alt.Scale(zero=False, padding=1)),
                            y=alt.Y('neg_log10_padj', title='-Log10 P-value', scale=alt.Scale(zero=False, padding=1)),
                            color=alt.Color('Regulation',
                                            scale=alt.Scale(domain=['Upregulated', 'Downregulated', 'Not Significant'],
                                                            range=['#ef4444', '#3b82f6', '#d1d5db'])),
                            tooltip=['log2FoldChange', 'padj', 'Regulation']
                        ).interactive().properties(height=500)
                        st.altair_chart(volcano_chart, use_container_width=True)
                    with c2:
                        st.markdown("**Top Upregulated**")
                        st.dataframe(
                            df_volcano[df_volcano['Regulation'] == 'Upregulated'].sort_values(by='padj').head(5)[
                                ['log2FoldChange', 'padj']], use_container_width=True)
                        st.markdown("**Top Downregulated**")
                        st.dataframe(
                            df_volcano[df_volcano['Regulation'] == 'Downregulated'].sort_values(by='padj').head(5)[
                                ['log2FoldChange', 'padj']], use_container_width=True)

    # ==========================================
    # VIEW 2: SINGLE SAMPLE DETAILS
    # ==========================================
    else:
        with st.sidebar:
            st.header("Drill Down")
            selected_sample = st.selectbox("Select Sample:", df_summary['Sample ID'].unique())

        st.markdown(f"### Deep Dive: `{selected_sample}`")
        tab_qc, tab_align, tab_exp = st.tabs(["🔍 FastP Detail", "🧬 Alignment Detail", "🧪 StringTie Detail"])

        with tab_qc:
            data_qc = load_fastp_detail(selected_sample)
            if data_qc:
                summ = data_qc['summary']['before_filtering']
                col1, col2, col3 = st.columns(3)
                col1.metric("Total Reads", f"{summ['total_reads'] / 1e6:.2f}M")
                col2.metric("Q30 Rate", f"{summ['q30_rate'] * 100:.2f}%")
                col3.metric("GC Content", f"{summ['gc_content'] * 100:.2f}%")
                col_q, col_c = st.columns(2)
                with col_q:
                    q_scores = data_qc['read1_before_filtering']['quality_curves']['mean']
                    df_q = pd.DataFrame({"Position": range(1, len(q_scores) + 1), "Score": q_scores}).set_index(
                        "Position")
                    chart_q = make_dynamic_line_chart(df_q, "Position", "Score", "Base Quality (Phred)", "Phred Score",
                                                      "#2563eb")
                    st.altair_chart(chart_q, use_container_width=True)
                with col_c:
                    curves = data_qc['read1_before_filtering']['content_curves']
                    df_c = pd.DataFrame(curves) * 100
                    df_c = df_c.reset_index().rename(columns={"index": "Position"})
                    df_c_melt = df_c.melt('Position', var_name='Base', value_name='Percentage')
                    chart_c = alt.Chart(df_c_melt).mark_line().encode(
                        x='Position', y=alt.Y('Percentage', scale=alt.Scale(zero=False, padding=1)), color='Base'
                    ).properties(title="Base Content (%)", height=300).interactive()
                    st.altair_chart(chart_c, use_container_width=True)

        with tab_align:
            h_metrics = load_hisat2_detail(selected_sample)
            if h_metrics:
                st.subheader("Alignment Statistics")
                c1, c2 = st.columns(2)
                c1.metric("Overall Alignment Rate", f"{h_metrics.get('alignment_rate', 0)}%")
                c1.metric("Total Reads Processed", f"{h_metrics.get('total_reads', 0):,}")
                df_align = pd.DataFrame([
                    {"Category": "Unique (1 time)", "Count": h_metrics.get('unique_aligned', 0)},
                    {"Category": "Multi (>1 time)", "Count": h_metrics.get('multi_aligned', 0)},
                    {"Category": "Unaligned", "Count": h_metrics.get('unaligned_concordant', 0)},
                ])
                with c2:
                    base = alt.Chart(df_align).encode(theta=alt.Theta("Count", stack=True))
                    pie = base.mark_arc(outerRadius=120, innerRadius=60).encode(
                        color=alt.Color("Category",
                                        scale=alt.Scale(domain=["Unique (1 time)", "Multi (>1 time)", "Unaligned"],
                                                        range=["#10b981", "#f59e0b", "#ef4444"])),
                        order=alt.Order("Count", sort="descending"), tooltip=["Category", "Count"]
                    )
                    st.altair_chart(pie, use_container_width=True)

        with tab_exp:
            df_st = load_stringtie_detail(selected_sample)
            if df_st is not None:
                st.subheader("🏆 Top 10 Most Expressed Genes")
                df_top = df_st.sort_values(by="TPM", ascending=False).head(10)
                st.table(df_top[['Gene ID', 'Gene Name', 'TPM', 'FPKM', 'Coverage']])
                st.bar_chart(df_top.set_index("Gene Name")["TPM"], color="#059669")
            else:
                st.warning("StringTie abundance file not found.")


if __name__ == "__main__":
    main()