"""Local research explorer for the qmsd result catalog."""
from __future__ import annotations
import json
import math

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from qmsd.asymptotics import optimal_gamma
from qmsd.puncture import column_to_point
from qmsd.results import (
    affine_line_profile, code_structure, distillation_series, load_result_catalog,
    pareto_front, puncture_points, record_dict,
)

STATUS_COLORS = {"confirmed": "#36d399", "partial": "#fbbf24",
                 "candidate": "#60a5fa", "refuted": "#fb7185"}
FIGURE_CONFIG = {"displaylogo": False, "toImageButtonOptions": {
    "format": "svg", "filename": "qmsd-figure", "scale": 2}}


def figure(fig):
    """Render with consistent publication-friendly export controls."""
    fig.update_layout(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(5,9,20,.55)", font={"family": "Inter, sans-serif"},
        margin={"l": 55, "r": 25, "t": 65, "b": 50})
    st.plotly_chart(fig, width="stretch", config=FIGURE_CONFIG)


def catalog_frame(records):
    return pd.DataFrame([{
        "artifact_id": r.artifact_id, "code": r.label, "status": r.status,
        "family": r.family, "p": r.p, "m": r.m, "n": r.n, "k": r.k,
        "d": r.d, "A_d": r.A_d, "gamma": r.gamma, "rate k/n": r.rate,
        "sublogarithmic": r.sublogarithmic, "provenance": r.provenance,
    } for r in records])


def metric(value):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "—"
    return f"{value:,}" if isinstance(value, int) else f"{value:.4g}"


def overview(records):
    st.subheader("Research landscape")
    confirmed = [r for r in records if r.status == "confirmed"]
    values = [len(records), len(confirmed), sum(r.sublogarithmic for r in confirmed),
              sum(r.A_d is not None for r in confirmed)]
    for col, title, value in zip(st.columns(4),
            ["Catalog artifacts", "Confirmed", "Sublogarithmic", "Exact A_d present"], values):
        col.metric(title, value)
    if not records:
        return
    frame = catalog_frame(records)
    left, right = st.columns([1.25, 1])
    with left:
        plot = frame.dropna(subset=["gamma"]).copy()
        if not plot.empty:
            plot["block size"] = plot.n.map(lambda x: max(8, min(38, 7 + math.log10(x + 1) * 7)))
            plot["Pareto"] = plot.artifact_id.isin(pareto_front(records, "gamma", "n"))
            fig = px.scatter(plot, x="n", y="gamma", color="status", symbol="Pareto",
                size="block size", hover_name="code", log_x=True,
                hover_data=["p", "m", "k", "d", "A_d", "family"],
                color_discrete_map=STATUS_COLORS, title="Finite-code yield landscape")
            fig.add_hline(y=1, line_dash="dash", line_color="#fbbf24",
                          annotation_text="gamma = 1")
            st.plotly_chart(fig, width="stretch")
    with right:
        counts = frame.groupby(["status", "family"], as_index=False).size()
        fig = px.bar(counts, x="size", y="family", color="status", orientation="h",
                     color_discrete_map=STATUS_COLORS, title="Evidence and construction mix")
        fig.update_layout(yaxis={"categoryorder": "total ascending"})
        figure(fig)

    primes = sorted({r.p for r in records})
    asymptotic = pd.DataFrame([{"p": p, "gamma": optimal_gamma(p)[0]} for p in primes])
    finite = frame[(frame.status == "confirmed") & frame.gamma.notna()]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=finite.p, y=finite.gamma, mode="markers",
        name="confirmed finite codes", text=finite.code,
        marker={"color": "#60a5fa", "size": 9, "opacity": .72}))
    fig.add_trace(go.Scatter(x=asymptotic.p, y=asymptotic.gamma, mode="lines+markers",
        name="asymptotic gamma₀(p)", line={"color": "#36d399", "width": 3}))
    fig.add_hline(y=1, line_dash="dash", line_color="#fbbf24")
    fig.update_layout(title="Finite results against the asymptotic yield limit",
                      xaxis_title="prime p", yaxis_title="gamma")
    st.plotly_chart(fig, width="stretch")


def structure_view(record):
    structure = code_structure(record)
    if structure is None:
        st.info("An explicit stabilizer/logical graph requires stored punctures and a tractable "
                "finite Reed–Muller block. This artifact is analytic-only or exceeds the display gate.")
        return

    x_stab, z_stab = structure["X_stab"], structure["Z_stab"]
    logical = structure["logical_X"]
    st.caption("Algebraic CSS incidence—not a hardware routing claim. A nonzero matrix entry "
               "connects a check or logical representative to a physical qudit.")
    for col, title, value in zip(st.columns(4),
            ["physical qudits", "X checks", "Z checks", "logical X representatives"],
            [record.n, x_stab.shape[0], z_stab.shape[0], logical.shape[0]]):
        col.metric(title, value)

    if logical.shape[0]:
        logical_index = st.selectbox("Logical qudit representative", range(logical.shape[0]),
                                     format_func=lambda i: f"logical X̄_{i + 1}",
                                     key=f"logical-{record.artifact_id}")
        coeff = logical[logical_index]
        coeff_by_column = dict(zip(structure["physical_columns"], coeff))
        ambient = []
        punctured = set(record.puncture_columns)
        for column in range(1, record.p ** record.m + 1):
            point = column_to_point(column, record.m, record.p)
            value = int(coeff_by_column.get(column, 0))
            role = "punctured" if column in punctured else (
                "logical support" if value else "physical")
            ambient.append((*point, column, value, role))
        dims = [f"x{i + 1}" for i in range(record.m)]
        data = pd.DataFrame(ambient, columns=dims + ["column", "coefficient", "role"])
        if record.m == 2:
            fig = px.scatter(data, x="x1", y="x2", color="role", symbol="coefficient",
                hover_data=["column", "coefficient"], title="Physical carrier of the selected logical X̄",
                color_discrete_map={"logical support": "#fbbf24", "physical": "#334155",
                                    "punctured": "#fb7185"})
        else:
            chosen = st.multiselect("Logical-map coordinates", dims,
                default=dims[:min(3, len(dims))], max_selections=3,
                key=f"logical-dims-{record.artifact_id}")
            if len(chosen) < 2:
                st.warning("Select at least two logical-map coordinates.")
                return
            args = dict(data_frame=data, color="role", symbol="coefficient",
                hover_data=["column", "coefficient"],
                color_discrete_map={"logical support": "#fbbf24", "physical": "#334155",
                                    "punctured": "#fb7185"})
            if len(chosen) == 2:
                fig = px.scatter(x=chosen[0], y=chosen[1], **args)
            else:
                fig = px.scatter_3d(x=chosen[0], y=chosen[1], z=chosen[2], **args)
            fig.update_layout(title="Physical carrier of the selected logical X̄")
        fig.update_traces(marker={"size": 7, "opacity": .82})
        figure(fig)

    st.markdown("#### Commutation by cancellation")
    st.caption("Each cell is the number of shared physical qudits between an X and Z check. "
               "The CSS condition is coefficient-weighted: the finite-field dot product "
               "vanishes even when the two supports overlap.")
    if x_stab.size and z_stab.size:
        support_overlap = (x_stab != 0).astype(int) @ (z_stab != 0).astype(int).T
        residual = x_stab @ z_stab.T % record.p
        max_x, max_z = min(35, support_overlap.shape[0]), min(35, support_overlap.shape[1])
        left, right = st.columns([1.35, 1])
        with left:
            overlap_fig = px.imshow(support_overlap[:max_x, :max_z], aspect="auto",
                origin="lower", color_continuous_scale="Magma",
                labels={"x": "Z-check", "y": "X-check", "color": "shared support"},
                title="X/Z support overlap (top-left check window)")
            figure(overlap_fig)
        with right:
            x_weights = np.count_nonzero(x_stab, axis=1)
            z_weights = np.count_nonzero(z_stab, axis=1)
            degrees = pd.DataFrame({"weight": np.r_[x_weights, z_weights],
                "check type": ["X"] * len(x_weights) + ["Z"] * len(z_weights)})
            degree_fig = px.histogram(degrees, x="weight", color="check type",
                barmode="overlay", opacity=.72,
                color_discrete_map={"X": "#60a5fa", "Z": "#fb7185"},
                title="Check-weight spectrum")
            figure(degree_fig)
        c1, c2, c3 = st.columns(3)
        c1.metric("largest X-check", int(x_weights.max(initial=0)))
        c2.metric("largest Z-check", int(z_weights.max(initial=0)))
        c3.metric("nonzero commutators", int(np.count_nonzero(residual)),
                  help="Must be zero for a valid CSS code.")

    matrix_name = st.radio("Incidence layer", ["X stabilizers", "Z stabilizers", "logical X"],
                           horizontal=True, key=f"matrix-{record.artifact_id}")
    matrix = {"X stabilizers": x_stab, "Z stabilizers": z_stab,
              "logical X": logical}[matrix_name]
    if matrix.size:
        max_rows = min(40, matrix.shape[0])
        shown = matrix[:max_rows]
        fig = px.imshow(shown, aspect="auto", color_continuous_scale="Turbo",
            labels={"x": "physical qudit index", "y": f"{matrix_name} row",
                    "color": f"F_{record.p} coefficient"},
            title=f"{matrix_name} incidence matrix"
                  + (f" (first {max_rows} rows)" if matrix.shape[0] > max_rows else ""))
        st.plotly_chart(fig, width="stretch")

        check_index = st.selectbox("Connectivity row", range(matrix.shape[0]),
            format_func=lambda i: f"{matrix_name} row {i + 1}",
            key=f"check-{record.artifact_id}-{matrix_name}")
        row = matrix[check_index]
        support = np.flatnonzero(row)
        if len(support):
            edge_x, edge_y = [], []
            ypos = np.linspace(-1, 1, len(support))
            for y in ypos:
                edge_x.extend([0, 1, None])
                edge_y.extend([0, y, None])
            graph = go.Figure()
            graph.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                line={"color": "#475569", "width": 1}, hoverinfo="skip"))
            graph.add_trace(go.Scatter(x=[0], y=[0], mode="markers+text",
                marker={"size": 28, "color": "#60a5fa"}, text=[f"{matrix_name} {check_index + 1}"],
                textposition="top center", name="row"))
            graph.add_trace(go.Scatter(x=np.ones(len(support)), y=ypos, mode="markers",
                marker={"size": 10, "color": row[support], "colorscale": "Turbo",
                        "cmin": 0, "cmax": record.p - 1, "showscale": True,
                        "colorbar": {"title": "coefficient"}},
                text=[f"physical {i + 1} · column {structure['physical_columns'][i]}"
                      for i in support], hovertemplate="%{text}<extra></extra>", name="physical"))
            graph.update_layout(title=f"Connectivity of one {matrix_name} row · weight {len(support)}",
                xaxis={"visible": False}, yaxis={"visible": False}, showlegend=False, height=430)
            figure(graph)

def geometry(record):
    pts = puncture_points(record)
    if not pts:
        st.info("No explicit puncture columns are registered for this artifact.")
        return
    data = pd.DataFrame(pts, columns=[f"x{i + 1}" for i in range(len(pts[0]))])
    data["column"] = list(record.puncture_columns)
    dims = list(data.columns[:record.m])
    data["Manhattan weight"] = data[dims].sum(axis=1)
    st.caption(f"{len(data):,} punctures in F_{record.p}^{record.m}; columns are 1-indexed.")
    if record.m == 2:
        fig = px.scatter(data, x="x1", y="x2", color="Manhattan weight",
                         hover_data=["column"], color_continuous_scale="Turbo",
                         title="Puncture geometry")
    elif record.m == 3:
        fig = px.scatter_3d(data, x="x1", y="x2", z="x3", color="Manhattan weight",
                            hover_data=["column"], color_continuous_scale="Turbo",
                            title="Puncture geometry")
    else:
        chosen = st.multiselect("Projection coordinates", dims,
            default=dims[:min(3, len(dims))], max_selections=3,
            key=f"dims-{record.artifact_id}")
        if len(chosen) < 2:
            st.warning("Select at least two coordinates.")
            return
        if len(chosen) == 2:
            fig = px.scatter(data, x=chosen[0], y=chosen[1], color="Manhattan weight",
                hover_data=["column"], color_continuous_scale="Turbo",
                title="Selected finite-field projection")
        else:
            fig = px.scatter_3d(data, x=chosen[0], y=chosen[1], z=chosen[2],
                color="Manhattan weight", hover_data=["column"],
                color_continuous_scale="Turbo", title="Selected finite-field projection")
    fig.update_traces(marker={"size": 7, "opacity": .82})
    figure(fig)
    hist = data[dims].melt(var_name="coordinate", value_name="value")
    counts = hist.groupby(["coordinate", "value"], as_index=False).size()
    heat = counts.pivot(index="coordinate", columns="value", values="size").fillna(0)
    figure(px.imshow(heat, text_auto=True, aspect="auto",
        color_continuous_scale="Viridis", title="Coordinate occupancy"))

    profile = affine_line_profile(record)
    if profile is not None:
        st.markdown("#### Affine-line obstruction scan")
        st.caption("Generalized Reed–Muller minimum-weight structure is affine-geometric. "
                   "This counts selected punctures on every affine line; a cap has at "
                   "most two on any line.")
        hist = pd.DataFrame(profile["histogram"],
                            columns=["punctures on line", "affine lines"])
        line_fig = px.bar(hist, x="punctures on line", y="affine lines", text_auto=True,
            color="punctures on line", color_continuous_scale="Turbo",
            title=f"Line-occupancy spectrum across {profile['line_count']:,} affine lines")
        line_fig.update_xaxes(dtick=1)
        figure(line_fig)
        max_occ = profile["max_occupancy"]
        verdict = "cap-like: no three collinear" if max_occ <= 2 else "contains collinear triples"
        st.metric("maximum punctures on one line", max_occ, verdict)
        witness = profile["max_line_columns"]
        st.caption("Example maximally occupied line (punctures highlighted): " + " · ".join(
            f"**{column}**" if column in set(record.puncture_columns) else str(column)
            for column in witness))


def distillation(record):
    low, high = st.slider("Input-error exponent range", 1, 6, (2, 5),
                          help="Plots delta_in from 10^-high to 10^-low.")
    df = pd.DataFrame(distillation_series(record, np.logspace(-high, -low, 100)))
    fig = go.Figure(go.Scatter(x=df.delta_in, y=df.cost, name="cost C",
                               line={"color": "#60a5fa"}))
    fig.update_xaxes(type="log", title="input error delta_in")
    fig.update_yaxes(type="log", title="inputs per accepted output")
    fig.update_layout(title="Single-round cost")
    figure(fig)
    if df.delta_out.notna().any():
        fig = px.line(df, x="delta_in", y="delta_out", log_x=True, log_y=True,
                      title="Leading-order output error")
        fig.update_traces(line={"color": "#36d399", "width": 3})
        figure(fig)
    else:
        st.info("Output-error curve requires confirmed distance and a known A_d.")


def dossier(record):
    color = "green" if record.status == "confirmed" else "orange"
    st.markdown(f"### {record.label} &nbsp; :{color}-badge[{record.status.upper()}]")
    st.caption(f"{record.family} · {record.provenance} · artifact {record.artifact_id}")
    values = [record.p, record.m, record.n, record.k, record.d, record.gamma, record.A_d]
    for col, name, value in zip(st.columns(7),
            ["p", "m", "n", "k", "d", "gamma", "A_d"], values):
        col.metric(name, metric(value))
    evidence, notes = st.columns([1.1, 1])
    with evidence:
        st.markdown("#### Certification ledger")
        st.write(f"**Distance:** {record.distance_evidence}")
        st.write(f"**A_d:** {record.Ad_evidence}")
        st.write(f"**Source:** {record.source}")
    with notes:
        st.markdown("#### Artifact note")
        st.write(record.note or "No additional note recorded.")
    tabs = st.tabs(["Code structure", "Puncture geometry", "Distillation", "Raw record"])
    with tabs[0]:
        structure_view(record)
    with tabs[1]:
        geometry(record)
    with tabs[2]:
        distillation(record)
    with tabs[3]:
        payload = record_dict(record)
        st.json(payload, expanded=False)
        st.download_button("Download normalized JSON", json.dumps(payload, indent=2),
            file_name=f"{record.artifact_id}.json", mime="application/json")
        if record.puncture_columns:
            st.download_button("Download puncture columns",
                "\n".join(map(str, record.puncture_columns)),
                file_name=f"{record.artifact_id}-punctures.txt")


def main():
    st.set_page_config(page_title="QMSD Research Explorer", page_icon="◈", layout="wide")
    st.markdown("""<style>
    .stApp { background:radial-gradient(circle at 80% 0%,#172554 0,#08111f 34%,#050914 75%); }
    [data-testid="stMetric"] { background:#0d1728;border:1px solid #20314d;padding:12px;border-radius:10px; }
    </style>""", unsafe_allow_html=True)
    st.title("QMSD Research Explorer")
    st.caption("Curated code properties, puncture geometry, distillation metrics, and evidence provenance")
    records = list(load_result_catalog())
    with st.sidebar:
        st.header("Catalog filters")
        statuses = st.multiselect("Evidence", list(STATUS_COLORS), default=["confirmed"])
        primes = st.multiselect("Prime p", sorted({r.p for r in records}))
        families = st.multiselect("Construction", sorted({r.family for r in records}))
        sublog = st.checkbox("gamma < 1 only")
    filtered = [r for r in records if r.status in statuses]
    if primes:
        filtered = [r for r in filtered if r.p in primes]
    if families:
        filtered = [r for r in filtered if r.family in families]
    if sublog:
        filtered = [r for r in filtered if r.sublogarithmic]
    overview(filtered)
    st.subheader("Code catalog")
    frame = catalog_frame(filtered)
    if not frame.empty:
        st.dataframe(frame.drop(columns=["artifact_id"]), width="stretch",
            hide_index=True, column_config={
                "gamma": st.column_config.NumberColumn(format="%.4f"),
                "rate k/n": st.column_config.NumberColumn(format="%.4f")})
        st.download_button("Export filtered catalog", frame.to_csv(index=False),
                           file_name="qmsd-results.csv", mime="text/csv")
    if not filtered:
        st.warning("No artifacts match the current filters.")
        return
    lookup = {f"{r.label} — {r.artifact_id}": r for r in filtered}
    compared = st.multiselect("Compare artifacts", list(lookup),
        default=list(lookup)[:min(3, len(lookup))], max_selections=6)
    if len(compared) >= 2:
        comparison = catalog_frame([lookup[name] for name in compared])
        comparison["cost @ 1e-3"] = [
            distillation_series(lookup[name], [1e-3])[0]["cost"] for name in compared
        ]
        st.markdown("#### Direct comparison")
        st.dataframe(comparison.drop(columns=["artifact_id", "provenance"]),
                     width="stretch", hide_index=True)
        chart = comparison.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["gamma", "cost @ 1e-3"])
        if not chart.empty:
            fig = px.scatter(chart, x="cost @ 1e-3", y="gamma", size="n", color="p",
                hover_name="code", hover_data=["k", "d", "A_d", "family"],
                log_x=True, color_continuous_scale="Turbo",
                title="Yield–cost comparison at delta_in = 10⁻³")
            fig.add_hline(y=1, line_dash="dash", line_color="#fbbf24")
            st.plotly_chart(fig, width="stretch")
    choices = list(lookup)
    structural_default = next((i for i, name in enumerate(choices)
        if lookup[name].puncture_columns and lookup[name].r is not None), 0)
    selected = st.selectbox("Inspect artifact", choices, index=structural_default, key="selected-artifact")
    dossier(lookup[selected])


if __name__ == "__main__":
    main()
