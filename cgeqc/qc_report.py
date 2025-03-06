import json
import os
from pathlib import Path
import base64
from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
from weasyprint import HTML, CSS
from jinja2 import Environment, FileSystemLoader
from scipy.stats import norm, lognorm

from cgeqc.qc_config import get_thresholds

def create_qc_report(trim_json_path, output_dir, name, pipeline_type="bacterial", trim_parameters=None):
    """Create a QC report from KMA trim output.
    
    Args:
        trim_json_path (str): Path to the trim json file
        output_dir (str): Output directory
        name (str): Run name
        pipeline_type (str): Type of data (bacterial, viral, or metagenomic)
        trim_parameters (dict, optional): Parameters used for trimming
    
    Returns:
        Path: Path to generated PDF report
    """
    # Load QC data
    with open(trim_json_path) as f:
        qc_data = json.load(f)
    
    # Calculate derived metrics
    metrics = calculate_qc_metrics(qc_data, pipeline_type)
    
    # Generate plots
    plots = generate_qc_plots(qc_data)
    
    # Create report
    return render_qc_report(metrics, plots, output_dir, name, trim_parameters)

def calculate_qc_metrics(qc_data, pipeline_type="bacterial"):
    """Calculate key QC metrics and determine quality assessment."""
    # Get appropriate thresholds
    thresholds = get_thresholds(pipeline_type)
    
    # Base metrics that are common for all types
    metrics = {
        'read_count': {
            'before': qc_data['Org. Fragment Count'],
            'after': qc_data['Fragment Count'],
            'change': calculate_percentage_change(
                qc_data['Fragment Count'], 
                qc_data['Org. Fragment Count']
            )
        },
        'total_bases': {
            'before': qc_data['Org. Bp Count'],
            'after': qc_data['Bp Count'],
            'change': calculate_percentage_change(
                qc_data['Bp Count'], 
                qc_data['Org. Bp Count']
            )
        },
        'mean_length': {
            'before': round(qc_data['Org. Mean Read Length'], 1),
            'after': round(qc_data['Mean Read Length'], 1),
            'change': calculate_percentage_change(
                qc_data['Mean Read Length'], 
                qc_data['Org. Mean Read Length']
            )
        },
        'n50': qc_data['N50'],
        'mean_quality': round(qc_data['E(Q)'], 1),
        'dataset_type': pipeline_type,
        'quality_assessment': {
            'status': 'good',
            'message': 'Overall data quality looks good',
            'points_to_check': []
        }
    }
    
    # Add type-specific metrics
    if pipeline_type == "bacterial":
        # Estimate coverage based on typical bacterial genome size
        TYPICAL_BACTERIAL_GENOME = 5_000_000  # 5 Mbp
        estimated_coverage = qc_data['Bp Count'] / TYPICAL_BACTERIAL_GENOME
        metrics['estimated_coverage'] = round(estimated_coverage, 1)
        metrics['gc_content'] = round(qc_data['GC Content'] * 100, 1)
    else:
        # For viral and metagenomic, we focus on total bp count
        metrics['bp_count'] = qc_data['Bp Count']
    
    # Quality checks and assessment
    assessment_points = []
    
    # Perform appropriate assessment based on data type
    if pipeline_type == "bacterial":
        # Bacterial assessment (coverage, quality, read length, GC content)
        quality_is_good = metrics['mean_quality'] >= thresholds['FAIR']['min_quality']
        coverage_is_good = metrics['estimated_coverage'] >= thresholds['FAIR']['min_coverage']
        quality_is_excellent = metrics['mean_quality'] >= thresholds['GOOD']['min_quality']
        coverage_is_excellent = metrics['estimated_coverage'] >= thresholds['GOOD']['min_coverage']
        
        # Base feedback on quality and coverage
        if not quality_is_good and not coverage_is_good:
            assessment_points.append(f'Both sequencing quality (Q{metrics["mean_quality"]:.1f}) and depth ({metrics["estimated_coverage"]:.1f}x) are below recommended levels. We recommend at least Q{thresholds["FAIR"]["min_quality"]} and {thresholds["FAIR"]["min_coverage"]}x coverage for reliable bacterial analysis')
        elif not quality_is_good:
            if coverage_is_excellent:
                assessment_points.append(f'Despite excellent sequencing depth ({metrics["estimated_coverage"]:.1f}x), the quality scores (Q{metrics["mean_quality"]:.1f}) are below recommended levels (Q{thresholds["FAIR"]["min_quality"]})')
            else:
                assessment_points.append(f'The quality scores (Q{metrics["mean_quality"]:.1f}) are below recommended levels (Q{thresholds["FAIR"]["min_quality"]}), though sequencing depth ({metrics["estimated_coverage"]:.1f}x) is adequate')
        elif not coverage_is_good:
            if quality_is_excellent:
                assessment_points.append(f'Despite excellent quality scores (Q{metrics["mean_quality"]:.1f}), the sequencing depth ({metrics["estimated_coverage"]:.1f}x) is below recommended levels ({thresholds["FAIR"]["min_coverage"]}x)')
            else:
                assessment_points.append(f'The sequencing depth ({metrics["estimated_coverage"]:.1f}x) is below recommended levels ({thresholds["FAIR"]["min_coverage"]}x), though quality scores (Q{metrics["mean_quality"]:.1f}) are adequate')
        else:
            # Both metrics are at least 'fair', give appropriate feedback
            if quality_is_excellent and coverage_is_excellent:
                assessment_points.append(f'Excellent sequencing quality (Q{metrics["mean_quality"]:.1f}) and depth ({metrics["estimated_coverage"]:.1f}x), both well above recommended levels')
            elif quality_is_excellent:
                assessment_points.append(f'Excellent quality scores (Q{metrics["mean_quality"]:.1f}). Sequencing depth ({metrics["estimated_coverage"]:.1f}x) is good but could be improved')
            elif coverage_is_excellent:
                assessment_points.append(f'Excellent sequencing depth ({metrics["estimated_coverage"]:.1f}x). Quality scores (Q{metrics["mean_quality"]:.1f}) are good but could be improved')
            else:
                assessment_points.append(f'Both quality scores (Q{metrics["mean_quality"]:.1f}) and sequencing depth ({metrics["estimated_coverage"]:.1f}x) are good, though could be improved for optimal results')
        
        # Read length checks
        if metrics['mean_length']['after'] < thresholds['POOR']['min_read_length']:
            assessment_points.append(f'The average read length is unusually short (below {thresholds["POOR"]["min_read_length"]} bp) for ONT bacterial sequencing')
        elif metrics['mean_length']['after'] < thresholds['FAIR']['min_read_length']:
            assessment_points.append(f'The average read length is shorter (< {thresholds["FAIR"]["min_read_length"]}) than typically seen with ONT bacterial sequencing')
        elif metrics['mean_length']['after'] >= thresholds['GOOD']['min_read_length']:
            assessment_points.append(f'Good average read length ({metrics["mean_length"]["after"]:.0f} bp) for ONT bacterial sequencing')
        
        # GC content checks for bacteria only
        if metrics['gc_content'] < thresholds['GC_CONTENT']['min'] or metrics['gc_content'] > thresholds['GC_CONTENT']['max']:
            assessment_points.append(f'The GC content ({metrics["gc_content"]}%) falls outside the typical range for bacterial genomes ({thresholds["GC_CONTENT"]["min"]}-{thresholds["GC_CONTENT"]["max"]}%). This might indicate potential contamination or bias in the sequencing')
        else:
            assessment_points.append(f'GC content ({metrics["gc_content"]}%) is within expected range for bacterial genomes ({thresholds["GC_CONTENT"]["min"]}-{thresholds["GC_CONTENT"]["max"]}%)')
        
        # Set overall assessment based on defined thresholds
        if (metrics['mean_quality'] >= thresholds['GOOD']['min_quality'] and 
            metrics['estimated_coverage'] >= thresholds['GOOD']['min_coverage']):
            metrics['quality_assessment'] = {
                'status': 'good',
                'message': 'Data quality is good for bacterial analysis',
                'points_to_check': assessment_points
            }
        elif (metrics['mean_quality'] >= thresholds['FAIR']['min_quality'] and 
              metrics['estimated_coverage'] >= thresholds['FAIR']['min_coverage']):
            metrics['quality_assessment'] = {
                'status': 'fair',
                'message': 'Data is suitable for bacterial analysis but some quality aspects may affect results',
                'points_to_check': assessment_points
            }
        else:
            metrics['quality_assessment'] = {
                'status': 'poor',
                'message': 'Data quality issues may significantly impact bacterial analysis reliability',
                'points_to_check': assessment_points
            }
    
    else:
        # Viral or metagenomic assessment (bp count and quality)
        bp_count = metrics['bp_count']
        quality_is_good = metrics['mean_quality'] >= thresholds['FAIR']['min_quality']
        bp_count_is_good = bp_count >= thresholds['FAIR']['min_bp_count']
        quality_is_excellent = metrics['mean_quality'] >= thresholds['GOOD']['min_quality']
        bp_count_is_excellent = bp_count >= thresholds['GOOD']['min_bp_count']
        
        data_type_name = "viral genome" if pipeline_type == "viral" else "metagenomic"
        
        # Base feedback on quality and bp_count
        if not quality_is_good and not bp_count_is_good:
            assessment_points.append(f'Both sequencing quality (Q{metrics["mean_quality"]:.1f}) and data volume ({bp_count:,} bp) are below recommended levels. We recommend at least Q{thresholds["FAIR"]["min_quality"]} and {thresholds["FAIR"]["min_bp_count"]:,} bp for reliable {data_type_name} analysis')
        elif not quality_is_good:
            if bp_count_is_excellent:
                assessment_points.append(f'Despite excellent data volume ({bp_count:,} bp), the quality scores (Q{metrics["mean_quality"]:.1f}) are below recommended levels (Q{thresholds["FAIR"]["min_quality"]})')
            else:
                assessment_points.append(f'The quality scores (Q{metrics["mean_quality"]:.1f}) are below recommended levels (Q{thresholds["FAIR"]["min_quality"]}), though data volume ({bp_count:,} bp) is adequate')
        elif not bp_count_is_good:
            if quality_is_excellent:
                assessment_points.append(f'Despite excellent quality scores (Q{metrics["mean_quality"]:.1f}), the data volume ({bp_count:,} bp) is below recommended levels ({thresholds["FAIR"]["min_bp_count"]:,} bp)')
            else:
                assessment_points.append(f'The data volume ({bp_count:,} bp) is below recommended levels ({thresholds["FAIR"]["min_bp_count"]:,} bp), though quality scores (Q{metrics["mean_quality"]:.1f}) are adequate')
        else:
            # Both metrics are at least 'fair', give appropriate feedback
            if quality_is_excellent and bp_count_is_excellent:
                assessment_points.append(f'Excellent sequencing quality (Q{metrics["mean_quality"]:.1f}) and data volume ({bp_count:,} bp), both well above recommended levels for {data_type_name} analysis')
            elif quality_is_excellent:
                assessment_points.append(f'Excellent quality scores (Q{metrics["mean_quality"]:.1f}). Data volume ({bp_count:,} bp) is good but could be improved for {data_type_name} analysis')
            elif bp_count_is_excellent:
                assessment_points.append(f'Excellent data volume ({bp_count:,} bp). Quality scores (Q{metrics["mean_quality"]:.1f}) are good but could be improved for {data_type_name} analysis')
            else:
                assessment_points.append(f'Both quality scores (Q{metrics["mean_quality"]:.1f}) and data volume ({bp_count:,} bp) are good, though could be improved for optimal {data_type_name} analysis')
        
        # Set overall assessment based on thresholds
        if (metrics['mean_quality'] >= thresholds['GOOD']['min_quality'] and 
            bp_count >= thresholds['GOOD']['min_bp_count']):
            metrics['quality_assessment'] = {
                'status': 'good',
                'message': f'Data quality is good for {data_type_name} analysis',
                'points_to_check': assessment_points
            }
        elif (metrics['mean_quality'] >= thresholds['FAIR']['min_quality'] and 
              bp_count >= thresholds['FAIR']['min_bp_count']):
            metrics['quality_assessment'] = {
                'status': 'fair',
                'message': f'Data is suitable for {data_type_name} analysis but some quality aspects may affect results',
                'points_to_check': assessment_points
            }
        else:
            metrics['quality_assessment'] = {
                'status': 'poor',
                'message': f'Data quality issues may significantly impact {data_type_name} analysis reliability',
                'points_to_check': assessment_points
            }
    
    return metrics

def calculate_percentage_change(new_value, old_value):
    """Calculate percentage change between two values."""
    if old_value == 0:
        return 0
    return abs(round(((new_value - old_value) / old_value) * 100, 1))

def generate_qc_plots(qc_data):
    """Generate QC plots with reference distributions."""
    plots = {}
    
    # Quality score distribution
    plt.figure(figsize=(10, 6))
    q_scores = range(len(qc_data['Q Distribution']))
    dist_data = qc_data['Q Distribution']
    
    # Find the last non-zero quality score
    max_q = max(i for i, v in enumerate(dist_data) if v > 0)
    
    # Normalize distribution for comparison
    total_reads = sum(dist_data[:max_q+1])
    norm_dist = [x/total_reads for x in dist_data[:max_q+1]]
    
    # Create reference distribution (normal distribution centered at Q15)
    reference_x = np.linspace(0, max_q, 100)
    reference_y = norm.pdf(reference_x, loc=15, scale=3)
    reference_y = reference_y / max(reference_y) * max(norm_dist)  # Scale to match data
    
    # Plot both distributions
    plt.bar(q_scores[:max_q+1], norm_dist, color='#4a90e2', alpha=0.6, label='Your Data')
    plt.plot(reference_x, reference_y, '--', color='#2ecc71', label='Typical Distribution', linewidth=2)
    
    plt.xlabel('Quality Score')
    plt.ylabel('Proportion of Reads')
    plt.title('Quality Score Distribution')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Set x-axis limits and ticks
    plt.xlim(-1, max_q + 1)
    tick_spacing = 5 if max_q > 30 else 2
    plt.xticks(range(0, max_q + 1, tick_spacing))
    
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    buf.seek(0)
    plots['quality_dist'] = base64.b64encode(buf.read()).decode('utf-8')
    
    # Read length distribution
    fig = plt.figure(figsize=(12, 8))
    
    # Create main plot
    ax_main = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    
    resolution = qc_data['Length Resolution']
    dist_data = np.array(qc_data['Length Distribution'])
    
    # Find the last non-zero bin
    max_bin = max(i for i, v in enumerate(dist_data) if v > 0)
    dist_data = dist_data[:max_bin + 1]  # Trim to relevant data
    
    # Calculate statistics on trimmed data
    cumsum = np.cumsum(dist_data)
    total_reads = cumsum[-1]
    percentile_99_idx = np.searchsorted(cumsum, 0.99 * total_reads)
    
    # Create bin edges and centers for plotting (only for the data we have)
    bin_edges = np.array([(i * resolution, (i + 1) * resolution) for i in range(len(dist_data))])
    bin_centers = bin_edges.mean(axis=1)
    
    # Plot main distribution (up to 99th percentile)
    ax_main.bar(bin_centers[:percentile_99_idx+1], 
                dist_data[:percentile_99_idx+1], 
                width=resolution * 0.9,
                color='#4a90e2', 
                alpha=0.6)
    
    # Calculate and plot statistics
    mean_length = qc_data['Mean Read Length']
    n50 = qc_data['N50']
    
    # Add vertical lines for mean and N50
    ax_main.axvline(mean_length, color='#2ecc71', linestyle='--', label=f'Mean ({mean_length:.0f} bp)')
    ax_main.axvline(n50, color='#e74c3c', linestyle='--', label=f'N50 ({n50:.0f} bp)')
    
    ax_main.set_xlabel('Read Length (bp)')
    ax_main.set_ylabel('Number of Reads')
    ax_main.set_title('Read Length Distribution (showing reads up to 99th percentile)')
    ax_main.grid(True, alpha=0.3)
    ax_main.legend()
    
    # Create reasonable tick spacing for main plot
    max_length_main = (percentile_99_idx + 1) * resolution
    desired_ticks = 6
    tick_size = max_length_main / desired_ticks
    magnitude = 10 ** np.floor(np.log10(tick_size))
    tick_size = np.round(tick_size / magnitude) * magnitude
    tick_positions = np.arange(0, max_length_main + tick_size, tick_size)
    ax_main.set_xticks(tick_positions)
    ax_main.set_xticklabels([f'{x/1000:.0f}k' if x >= 1000 else str(int(x)) for x in tick_positions])
    
    # Create small overview plot
    ax_overview = plt.subplot2grid((3, 3), (2, 0), colspan=3)
    ax_overview.bar(bin_centers, dist_data, width=resolution * 0.9, color='#4a90e2', alpha=0.6)
    ax_overview.axvline(mean_length, color='#2ecc71', linestyle='--')
    ax_overview.axvline(n50, color='#e74c3c', linestyle='--')
    ax_overview.set_xlabel('Read Length (bp)')
    ax_overview.set_ylabel('Reads')
    ax_overview.set_title('Full Range Distribution')
    
    # Set x-axis for overview plot
    max_length = len(dist_data) * resolution
    tick_positions = np.linspace(0, max_length, 6)
    ax_overview.set_xticks(tick_positions)
    ax_overview.set_xticklabels([f'{x/1000:.0f}k' if x >= 1000 else str(int(x)) for x in tick_positions])
    
    # Add key statistics
    stats_text = (f'Mean: {mean_length:.0f} bp\n'
                 f'N50: {n50:.0f} bp\n'
                 f'Bin size: {resolution} bp\n'
                 f'Total reads: {int(total_reads):,}')
    
    ax_main.text(0.98, 0.98, stats_text,
                transform=ax_main.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    buf.seek(0)
    plots['length_dist'] = base64.b64encode(buf.read()).decode('utf-8')
    
    return plots

def render_qc_report(metrics, plots, output_dir, name, trim_parameters=None):
    """Render the QC report using the template."""
    package_dir = Path(__file__).parent
    env = Environment(
        loader=FileSystemLoader(package_dir / "templates"),
        autoescape=True
    )
    
    # Load logo
    logo_path = package_dir / "assets" / "dtu_logo.png"
    logo_data = ""
    if logo_path.exists():
        with open(logo_path, "rb") as f:
            logo_data = f'data:image/png;base64,{base64.b64encode(f.read()).decode("utf-8")}'
    
    # Ensure parameters are integers for correct template comparison
    if trim_parameters:
        for key in trim_parameters:
            if isinstance(trim_parameters[key], (int, float)):
                trim_parameters[key] = int(trim_parameters[key])
    
    # Get the pipeline type directly from metrics
    pipeline_type = metrics.get('dataset_type', 'bacterial')
    
    template = env.get_template("qc_report.html")
    html_content = template.render(
        name=name,
        metrics=metrics,
        plots=plots,
        logo_data_url=logo_data,
        pipeline_type=pipeline_type,
        trim_parameters=trim_parameters or {}  # Provide default empty dict if None
    )
    
    # Create PDF
    pdf_path = Path(output_dir) / f"{name}_qc_report.pdf"
    HTML(string=html_content).write_pdf(
        pdf_path,
        stylesheets=[CSS(package_dir / "assets" / "style.css")]
    )
    
    return pdf_path

