var genome_coords_margin = {left: 80, right: 30};
var coverage_plot_height = 100;
var coverage_plot_margin = {top: 10, bottom: 20};
var gene_plot_height = 15;
var gene_plot_margin = {top: 0, bottom: 0};
var variant_plot_height = 15;
var variant_plot_margin = {top: 0, bottom: 0};

function bootstrap_plot() {
    window.model = window.model || {};
    window.model.plot = window.model.plot || {};
    if (!window.model.plot.svg_width) {
        window.model.plot.svg_width = $('#coverage_plot_container').width();
        window.model.plot.genome_coords_width = window.model.plot.svg_width - genome_coords_margin.left - genome_coords_margin.right;
    }
    if (!window.model.plot.x) {
        window.model.plot.x = d3.scale.linear()
            .domain([window.model.start, window.model.stop])
            .range([0, window.model.plot.genome_coords_width]);
    }
}

function create_coverage_plot() {
    var coverage_XHR = $.getJSON(fmt('/api/coverage/region/{0}-{1}-{2}', window.model.chrom, window.model.start, window.model.stop));

    $(function() {
        bootstrap_plot();

        var xAxis = d3.svg.axis()
            .scale(window.model.plot.x)
            .orient("bottom")
            .ticks(5)
            .tickFormat(group_thousands);

        var svg = d3.select('#coverage_plot_container').append("svg")
            .attr("width", window.model.plot.svg_width)
            .attr("height", coverage_plot_height + coverage_plot_margin.top + coverage_plot_margin.bottom);
        var genome_g = svg.append("g")
            .attr('id', 'inner_graph')
            .attr('class', 'genome_g')
            .attr("transform", "translate(" + genome_coords_margin.left + "," + coverage_plot_margin.top + ")");

        genome_g.append('clipPath')
            .attr('id', 'cov-plot-clip')
            .append('rect')
            .attr('x', 0)
            .attr('width', window.model.plot.genome_coords_width)
            .attr('y', 0)
            .attr('height', coverage_plot_height);

        genome_g.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + coverage_plot_height + ")")
            .call(xAxis);

        var cov_bar_g = genome_g.append('g')
            .attr('clip-path', 'url(#cov-plot-clip)')
            .attr('id', 'cov-bar-g');

        genome_g.append('rect').attr('class', 'mouse_guide')
            .attr('x', -999).style('height', '100%').style('width', '1px')
            .style('fill', 'rgb(210,210,210)').style('fill-opacity', '85%');

        var loading_text = genome_g.append('text')
            .attr('text-anchor','middle').text('loading...')
            .attr('transform', fmt('translate({0},{1})', window.model.plot.genome_coords_width/2, coverage_plot_height/2))

        coverage_XHR
            .done(function(coverage_stats) {
                loading_text.remove();
                window.model.coverage_stats = coverage_stats;
                if (window.model.coverage_stats !== null) populate_coverage_plot(cov_bar_g, genome_g);
            })
            .fail(function() { console.error('coverage XHR failed'); });
    });
}

function populate_coverage_plot(cov_bar_g, genome_g) {
    var metric = 'mean';
    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(window.model.coverage_stats, function(d) { return d[metric]; });
    }
    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([coverage_plot_height, 0]);
    var yAxis = _coverage_y_axis(y, metric);

    cov_bar_g
        .selectAll("rect.cov_plot_bars")
        .data(window.model.coverage_stats)
        .enter()
        .append("rect")
        .attr('class', 'cov_plot_bars')
        .style("fill", "steelblue")
        .attr("x", function(d) {
            return window.model.plot.x(d.start);
        })
        .attr("width", function(d) {
            return window.model.plot.x(d.end) - window.model.plot.x(d.start) + 1;
        })
        .attr("y", function(d) {
            return y(d[metric]) || 0;
        })
        .attr("height", function(d) {
            return (coverage_plot_height - y(d[metric])) || 0;
        });

    genome_g.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    // Handle changes
    $('.coverage_metric_buttons').change(function () {
        var v = $(this).attr('id').replace('_covmet_button', '');
        $('.coverage_subcat_selectors').hide();
        if (v == 'covered') {
            $('#over_x_select_container').show();
            v = $('#over_x_select').val().replace('X', '');
        } else {
            $('#average_select_container').show();
            v = $("#average_select").val();
        }
        change_coverage_plot_metric(v);
    });
    $('#over_x_select').change(function () {
        change_coverage_plot_metric($(this).val().replace('X', ''));
    });
    $('#average_select').change(function () {
        change_coverage_plot_metric($(this).val());
    });
}

function _coverage_y_axis(y_scale, metric) {
    var yAxis = d3.svg.axis()
        .scale(y_scale)
        .orient('left')
        .ticks(3);

    if (metric === 'mean' || metric === 'median')
        yAxis = yAxis.tickFormat(function(d) {return d.toString() + '\u00d7'});
    else
        yAxis = yAxis.tickFormat(d3.format('%'));

    return yAxis;
}

function change_coverage_plot_metric(metric) {
    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(window.model.coverage_stats, function(d) { return d[metric]; });
    }

    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([coverage_plot_height, 0]);

    var svg = d3.select('#coverage_plot_container').select('#inner_graph');

    var yAxis = _coverage_y_axis(y, metric);
    svg.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);

    svg.select('#cov-bar-g')
        .selectAll("rect.cov_plot_bars")
        .data(window.model.coverage_stats)
        .transition()
        .duration(500)
        .attr("y", function(d) { return y(d[metric]); })
        .attr("height", function(d) { return coverage_plot_height - y(d[metric]); });
}

function create_gene_plot() {
    bootstrap_plot();

    var svg = d3.select('#gene_plot_container').append('svg')
        .attr('width', window.model.plot.svg_width)
        .attr('height', gene_plot_height)
    var genome_g = svg.append('g')
        .attr('class', 'genome_g')
        .attr('id', 'gene_track')
        .attr('transform', 'translate(' + genome_coords_margin.left+','+0+')');

    genome_g.append('rect').attr('class', 'mouse_guide')
        .attr('x', -999).style('height', '100%').style('width', '1px')
        .style('fill', 'rgb(210,210,210)').style('fill-opacity', '85%');

    window.model.plot.exon_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        return (d.feature_type==='CDS'?'Coding Sequence':'UTR') + '<br>' +
            'start: ' + group_thousands_html(d.start) + '<br>' +
            'stop: ' + group_thousands_html(d.stop);
    });
    svg.call(window.model.plot.exon_tip);

    genome_g.append('line')
        .attr("y1", gene_plot_height/2)
        .attr("y2", gene_plot_height/2)
        .attr("x1", 0)
        .attr("x2", window.model.plot.genome_coords_width)
        .attr("stroke-width", 1)
        .attr("stroke", "lightsteelblue");

    genome_g.selectAll('rect.exon')
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
}

function create_variant_plot() {
    bootstrap_plot();

    var svg = d3.select('#variant_plot_container').append("svg")
        .attr("width", window.model.plot.svg_width)
        .attr("height", variant_plot_height);
    var genome_g = svg.append("g")
        .attr('id', 'variant_track')
        .attr('class', 'genome_g')
        .attr("transform", "translate(" + genome_coords_margin.left + "," + 0 + ")");

    genome_g.append('rect').attr('class', 'mouse_guide')
        .attr('x', -999).style('height', '100%').style('width', '1px')
        .style('fill', 'rgb(210,210,210)').style('fill-opacity', '85%');

    genome_g.append("line")
        .attr("y1", variant_plot_height/2)
        .attr("y2", variant_plot_height/2)
        .attr("x1", 0)
        .attr("x2", window.model.plot.genome_coords_width)
        .attr("stroke-width", 1)
        .attr("stroke", "lightsteelblue")
        .style('opacity', 0.3);

    genome_g.append("line")
        .attr('id', 'variant_plot_region_selector')
        .attr("y1", variant_plot_height/2)
        .attr("y2", variant_plot_height/2)
        .attr("stroke-width", 1)
        .attr("stroke", "black");
    change_variant_plot_region_selector(window.model.start, window.model.stop)

    window.model.plot.oval_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        var csq = d.major_consequence.replace(/_/g, ' ');
        if (csq.length > 15) { csq = csq.substr(0, 15) + '...'; } // because d3-tip tooltips fall off the page
        return group_thousands_html(d.pos) + '<br>' +
            csq + '<br>' +
            (d.filter === 'PASS' ? '' : 'FAIL<br>') +
            'MAF: ' + perc_sigfigs_html(d.allele_freq, 2);
        //return JSON.stringify(d);
    });
    svg.call(window.model.plot.oval_tip);
}

function change_variant_plot_region_selector(start, stop) {
    var x1 = window.model.plot.x(start);
    var x2 = window.model.plot.x(stop);
    d3.select('#variant_plot_region_selector')
        .attr('x1', x1)
        .attr('x2', Math.max(x1, x2));
}

function change_variant_plot(variants) {
    var selection = d3.select('#variant_track')
        .selectAll('.variant-circle')
        .data(variants, get_variant_id); // define data-joining 2nd method to allow move-to-front
    selection.enter()
        .append('ellipse')
        .attr('class', 'variant-circle')
        .style('opacity', 0.3)
        .style('fill', 'blue')
        .attr('ry', 6)
        .attr('rx', 6)
        .attr('cy', variant_plot_height/2)
        .attr('cx', function(d) { return window.model.plot.x(d.pos); })
        .attr('id', get_variant_plot_id)
        .on('mouseover', function(variant) {
            window.model.plot.oval_tip.show(variant);
            $('.variant-circle').css('fill', 'blue').css('opacity', 0.3);
            $(this).css('fill', 'orange').css('opacity', 1);
            window.model.tbl.rows()[0].forEach(function(row_idx) {
                $(window.model.tbl.row(row_idx).nodes()).removeClass('highlight');
            });
            window.model.tbl.rows()[0].forEach(function(row_idx) {
                var cur_var = window.model.tbl.row(row_idx).data();
                if (cur_var.pos === variant.pos && cur_var.ref === variant.ref && cur_var.alt === variant.alt) {
                    $(window.model.tbl.row(row_idx).nodes()).addClass('highlight');
                }
            });
            d3.select(this).moveToFront();
        })
        .on('mouseout', function(variant) {
            $(this).css('fill', 'blue').css('opacity', 0.3);
            window.model.plot.oval_tip.hide(variant);
            window.model.tbl.rows()[0].forEach(function(row_idx) {
                var cur_var = window.model.tbl.row(row_idx).data();
                if (cur_var.pos === variant.pos && cur_var.ref === variant.ref && cur_var.alt === variant.alt) {
                    $(window.model.tbl.row(row_idx).nodes()).removeClass('highlight');
                }
            });
        })
    selection.exit()
        .remove()
}


function create_variant_table() {
    var columns = [
        {
            title: 'Alleles (rsID)', name: 'allele',
            searchable: false, orderable: false,
            render: (function() {
                var template = _.template(
                    '<a href="/variant/<%= window.model.chrom %>-<%= variant.pos %>-<%= variant.ref %>-<%= variant.alt %>" target="_blank">'+
                        '<% if (variant.ref.length > 20) { %><%= variant.ref.slice(0,20) %>...<% } else { %><%= variant.ref %><% } %> / '+
                        '<% if (variant.alt.length > 20) { %><%= variant.alt.slice(0,20) %>...<% } else { %><%= variant.alt %><% } %>'+
                        '<% if (variant.rsids.length) { %> (<%= variant.rsids.join(", ") %>)<% } %>'+
                        '</a>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Position (chr' + window.model.chrom + ')', name: 'pos',
            data: 'pos', searchable: true,  orderable: true, className: 'dt-right',
            render: function(cell_data, type, row) { return group_thousands_html(cell_data); },

        },{
            title: 'HGVS', name: 'hgvs',
            searchable: false, orderable: false, className: 'dt-head-center',
            render: (function() {
                var template = _.template(
                    '<% if (variant.HGVS.length > 20) { %><%= variant.HGVS.slice(0,20) %>...<% } else { %><%= variant.HGVS %><% } %>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Annotation', name: 'csq',
            data: 'major_consequence', searchable:true, orderable:false, className: 'dt-center',
            render: function(cell_data, type, row) { return fmt_annotation(cell_data); },

        },{
            title: 'CADD', name:'cadd_phred',
            data: 'cadd_phred', orderable:true, orderSequence:['desc', 'asc'], className: 'dt-right',
            render: function(cell_data, type, row) { return (cell_data===null)?'':cell_data.toFixed(0); },

        },{
            title: 'QC', name: 'filter',
            data: 'filter', searchable:true, orderable:false, className: 'dt-center',
            render: function(cell_data, type, row) { return (cell_data==='PASS') ? 'PASS' : fmt('<span data-tooltip="failed filters: {0}">FAIL</span>', cell_data); },

        },{
            title: 'N Alleles', name: 'allele_num',
            data: 'allele_num', searchable:true, orderable:true, className: 'dt-right',
            render: function(cell_data, type, row) {return group_thousands_html(cell_data);},

        },{
            title: 'Het', name: 'het',
            data: 'het', searchable:true, orderable:true, orderSequence:['desc','asc'], className: 'dt-right',
            render: function(cell_data, type, row) {return group_thousands_html(cell_data);},

        },{
            title: 'HomAlt', name: 'hom_count',
            searchable:true, orderable:true, orderSequence:['desc','asc'], className: 'dt-right',
            render: function(cell_data, type, row) {return group_thousands_html(row.hom_count);},

        },{
            title: 'Frequency', name: 'allele_freq',
            data: 'allele_freq', searchable:true, orderable:true, orderSequence:['desc','asc'], className: 'dt-pad-left',
            render: function(cell_data, type, row) { return perc_sigfigs_html(cell_data, 2); },

        },
    ];

    window.model = window.model || {};
    window.model.filter_info = window.model.filter_info || {};
    window.model.filter_info.start = window.model.start;
    window.model.filter_info.stop = window.model.stop;
    window.model.filter_info.chrom = window.model.chrom;

    window.model.tbl = $('#variant_table').DataTable({
        serverSide: true, /* API does all the real work */

        processing: true, /* show "processing" over table while waiting for API */
        deferRender: true, /* only render rows when they're being displayed */

        paging: true,
        pagingType: 'full', /* [first, prev, next, last] */
        pageLength: 10,

        searching: false,

        ajax: {
            url: '/api/variants_for_table',
            type: 'POST',
            data: function(args) { /* modify request form parameters */
                return {
                    args: JSON.stringify(args), // jsonify all params rather than using `columns[0][search][value]` php form syntax
                    filter_info: JSON.stringify(window.model.filter_info),
                };
            },
            dataSrc: function(resp) { /* modify API's response (for debugging) */
                window._debug = window._debug || {}; window._debug.resp = resp;
                change_variant_plot(resp.data);
                return resp.data;
            },
        },

        order: [[columns.map(function(d){return d.name==='cadd_phred'}).indexOf(true), 'desc']],
        columns: columns,

        dom: '<ipl>rft', // default is 'lfrtip'.  l=length f=filtering t=table i=info p=paging, r=processing

        language: {
            info: 'Showing variants _START_ - _END_ of _TOTAL_',
            infoFiltered: '(filtered from _MAX_ variants)',
            infoEmpty: 'No matching variants',
            thousands: '\u202f',
            lengthMenu: 'Show _MENU_ variants',
        }

    });

    // hilite corresponding variant-plot circle
    $('#variant_table tbody').on('mouseleave', 'tr', function() {
        var variant = window.model.tbl.row(this).data();
        if (variant) {
            $('#' + get_variant_plot_id(variant))
            .css('fill', 'blue')
            .css('opacity', 0.3);
        }
    });
    $('#variant_table tbody').on('mouseenter', 'tr', function() {
        var variant = window.model.tbl.row(this).data();
        $('.variant-circle')
            .css('fill', 'blue')
            .css('opacity', 0.3);
        if (variant) {
            var vid = '#' + get_variant_plot_id(variant);
            d3.select(vid).moveToFront();
            $(vid)
                .css('fill', 'orange')
                .css('opacity', 1);
        }
    });

    $('.variant_table_filter').on('change', function() {
        window.model.filter_info.pos_ge = parseInt($('input#pos_ge').val());
        window.model.filter_info.pos_le = parseInt($('input#pos_le').val());
        window.model.filter_info.maf_ge = parseFloat($('input#maf_ge').val()) / 100; // %
        window.model.filter_info.maf_le = parseFloat($('input#maf_le').val()) / 100; // %
        window.model.filter_info.filter_value = $('select#filter_value').val();
        window.model.tbl.draw();
    });

    $('input#pos_le,input#pos_ge').change(function() {
        change_variant_plot_region_selector(parseInt($('input#pos_ge').val()), parseInt($('input#pos_le').val()));
    });
}

create_coverage_plot();
$(function() {
    $('#pos_ge').val(window.model.start);
    $('#pos_le').val(window.model.stop);
    $('#pos_ge,#pos_le').attr('min', window.model.start);
    $('#pos_ge,#pos_le').attr('max', window.model.stop);
    $('#pos_ge,#pos_le').attr('step', Math.ceil((window.model.stop - window.model.start) / 20));
    create_gene_plot();
    create_variant_plot();
    create_variant_table();
    d3.selectAll('.genome_g')
        .on('mousemove', function() {
            var coords = d3.mouse(d3.select(this).node());
            d3.selectAll('.mouse_guide').attr('x', coords[0]);
        })
        .on('mouseleave', function() {
            d3.selectAll('.mouse_guide').attr('x', -999);
        })
        .append('rect').style('width', '100%').style('height', '100%').style('opacity', 0); // recives mousemove and bubbles it up to `.genome_g`
});

function get_variant_plot_id(variant) { return 'variant-plot-'+get_variant_id(variant); }
function get_variant_id(variant) { return ''+variant.pos+'-'+variant.ref+'-'+variant.alt; }
