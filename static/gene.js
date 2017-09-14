var gene_plot_height = 15;
var gene_plot_margin = {top: 0, bottom: 0};

function create_gene_plot() {
    bootstrap_plot();

    var svg = d3.select('#gene_plot_container').append('svg')
        .attr('width', window.model.plot.svg_width)
        .attr('height', gene_plot_height)
    var genome_g = svg.append('g')
        .attr('class', 'genome_g')
        .attr('id', 'gene_track')
        .attr('transform', 'translate(' + genome_coords_margin.left+','+0+')');

    genome_g.append('clipPath')
        .attr('id', 'gene-plot-clip')
        .append('rect')
        .attr('x', 0)
        .attr('width', window.model.plot.genome_coords_width)
        .attr('y', 0)
        .attr('height', gene_plot_height);

    window.model.plot.exon_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        return (d.feature_type==='CDS'?'Coding Sequence':'UTR') + '<br>' +
            'start: ' + group_thousands_html(d.start) + '<br>' +
            'stop: ' + group_thousands_html(d.stop);
    });
    svg.call(window.model.plot.exon_tip);

    var data_g = genome_g.append('g');

    data_g.selectAll('line.intervals')
        .data(window.model.intervalset.list_of_pairs)
        .enter()
        .append('line')
        .attr("y1", gene_plot_height/2)
        .attr("y2", gene_plot_height/2)
        .attr("x1", function(d) { return window.model.plot.x(d[0]) })
        .attr("x2", function(d) { return window.model.plot.x(d[1]) })
        .attr("stroke-width", 1)
        .attr("stroke", "lightsteelblue");

    data_g.selectAll('rect.exon')
        .data(window.model.exons)
        .enter()
        .append('rect')
        .attr('class', 'exon')
        .style('fill', 'lightsteelblue')
        .attr('y', function(d){return d.feature_type==='CDS' ? 0 : gene_plot_height/4})
        .attr('height',function(d){return d.feature_type==='CDS' ? gene_plot_height : gene_plot_height/2})
        .attr('x', function(d) { return window.model.plot.x(d.start) })
        .attr('width', function(d) { return window.model.plot.x(d.stop)-window.model.plot.x(d.start) })
        .on('mouseover', window.model.plot.exon_tip.show)
        .on('mouseout', window.model.plot.exon_tip.hide)

    genome_g.append('rect').attr('class', 'mouse_guide').attr('x', -999).attr('clip-path', 'url(#gene-plot-clip)');
}

create_summary_table();
create_coverage_plot();
$(function() {
    create_gene_plot();
    create_variant_plot();
    create_variant_table();
    initiate_mouse_guide();
});
