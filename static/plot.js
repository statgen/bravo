var genome_coords_margin = {left: 80, right: 30};

var coverage_plot_height = 70;
var coverage_plot_margin = {top: 8, bottom: 10};

var gene_plot_height = 15;
var gene_plot_margin = {top: 0, bottom: 3};

var pos_plot_height = 15;
var pos_plot_margin = {top: 0, bottom: 0};

var variant_plot_height = 20;
var variant_plot_margin = {top: 0, bottom: 0};


function bootstrap_plot() {
    window.model = window.model || {};
    window.model.plot = window.model.plot || {};
    if (!window.model.plot.svg_width) {
        window.model.plot.svg_width = $('#coverage_plot_container').width();
        window.model.plot.genome_coords_width = window.model.plot.svg_width - genome_coords_margin.left - genome_coords_margin.right;
    }
    if (!window.model.plot.x) {
        var total_range_length = sum(window.model.intervalset.list_of_pairs.map(function(pair){return pair[1]-pair[0]}));
        var SPACING = total_range_length/50;
        var domain = [], range = [];
        window.model.intervalset.list_of_pairs.forEach(function(pair) {
            var range_start = (range.length===0 ? 0 : range[range.length-1]+SPACING);
            var range_end = range_start + pair[1]-pair[0];
            domain.push(pair[0]); range.push(range_start);
            domain.push(pair[1]); range.push(range_end);
        });
        range = range.map(function(x) { return x * window.model.plot.genome_coords_width / range[range.length-1]; });
        window.model.plot.x = d3.scale.linear().domain(domain).range(range);
    }
}

function get_variant_plot_id(variant) { return 'variant-plot-'+get_variant_id(variant); }
function get_variant_id(variant) { return ''+variant.pos+'-'+variant.ref+'-'+variant.alt; }


function create_coverage_plot() {
    var coverage_XHR = $.getJSON(window.model.url_prefix + 'api/coverage' + window.model.url_suffix);
    $(function() {
        bootstrap_plot();

        var svg = d3.select('#coverage_plot_container').append("svg")
            .attr("width", window.model.plot.svg_width)
            .attr("height", coverage_plot_height + coverage_plot_margin.top + coverage_plot_margin.bottom)
            .style('display', 'block');
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

        var data_g = genome_g.append('g')
            .attr('clip-path', 'url(#cov-plot-clip)')
            .attr('id', 'cov-bar-g');

        genome_g.append('rect').attr('class', 'mouse_guide').attr('x', -999).attr('clip-path', 'url(#cov-plot-clip)')

        var loading_text = genome_g.append('text')
            .attr('text-anchor','middle').text('loading...')
            .attr('transform', fmt('translate({0},{1})', window.model.plot.genome_coords_width/2, coverage_plot_height/2))

        coverage_XHR
            .done(function(coverage_stats) {
                loading_text.remove();
                window.model.coverage_stats = coverage_stats;
                if (window.model.coverage_stats !== null) populate_coverage_plot(data_g, genome_g);
            })
            .fail(function() { console.error('coverage XHR failed'); });
    });
}

function populate_coverage_plot(data_g, genome_g) {
    var metric = 'mean';
    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(window.model.coverage_stats, function(d) { return d[metric]; });
    }
    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([coverage_plot_height, 0]);
    var yAxis = _coverage_y_axis(y, metric);

    data_g
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
            return window.model.plot.x(d.end + 1) - window.model.plot.x(d.start) + 1;
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

    window.model.exons = [].concat( // put UTR before CDS, so that UTR is rendered behind CDS
        window.model.exons.filter(function(exon) { return exon.feature_type !== 'CDS' }),
        window.model.exons.filter(function(exon) { return exon.feature_type === 'CDS' })
    );

    var svg = d3.select('#gene_plot_container').append('svg')
        .attr('width', window.model.plot.svg_width)
        .attr('height', gene_plot_height + gene_plot_margin.top + gene_plot_margin.bottom)
        .style('display', 'block');
    var genome_g = svg.append('g')
        .attr('class', 'genome_g')
        .attr('id', 'gene_track')
        .attr('transform', fmt('translate({0},{1})', genome_coords_margin.left, gene_plot_margin.top));

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


function create_pos_plot() {
    bootstrap_plot();

    var svg = d3.select('#pos_plot_container').append('svg')
        .attr('width', window.model.plot.svg_width)
        .attr('height', pos_plot_height + pos_plot_margin.top + pos_plot_margin.bottom)
        .style('display', 'block');
    var genome_g = svg.append('g')
        .attr('class', 'genome_g')
        .attr('transform', fmt('translate({0},{1})', genome_coords_margin.left, pos_plot_margin.top));
    genome_g.append('text').attr('id', 'pos_plot_text')
        .attr('text-anchor', 'middle');
}



function create_variant_plot() {
    bootstrap_plot();

    var svg = d3.select('#variant_plot_container').append("svg")
        .attr("width", window.model.plot.svg_width)
        .attr("height", variant_plot_height + variant_plot_margin.top + variant_plot_margin.bottom)
        .style('display', 'block');
    var genome_g = svg.append("g")
        .attr('id', 'variant_track')
        .attr('class', 'genome_g')
        .attr("transform", fmt('translate({0},{1})', genome_coords_margin.left, variant_plot_margin.top));

    genome_g.append('clipPath')
        .attr('id', 'variant-plot-clip')
        .append('rect')
        .attr('x', 0)
        .attr('width', window.model.plot.genome_coords_width)
        .attr('y', 0)
        .attr('height', coverage_plot_height);

    genome_g.append('rect').attr('class', 'mouse_guide').attr('x', -999).attr('clip-path', 'url(#variant-plot-clip)');

    genome_g.selectAll('line.intervals')
        .data(window.model.intervalset.list_of_pairs)
        .enter()
        .append('line')
        .attr("y1", variant_plot_height/2)
        .attr("y2", variant_plot_height/2)
        .attr("x1", function(d) { return window.model.plot.x(d[0]) })
        .attr("x2", function(d) { return window.model.plot.x(d[1]) })
        .attr("stroke-width", 1)
        .attr("stroke", "lightsteelblue");

    var data_g = genome_g.append('g').attr('class','data_g');

    window.model.plot.oval_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        return group_thousands_html(d.pos) + '<br>' +
            fmt_annotation(d, 20) + '<br>' +
            (d.filter === 'PASS' ? '' : 'FAIL<br>') +
            'MAF: ' + perc_sigfigs_html(d.allele_freq, 2);
        //return JSON.stringify(d);
    });
    svg.call(window.model.plot.oval_tip);
}

function change_variant_plot(variants) {
    var selection = d3.select('#variant_track .data_g')
        .selectAll('.variant-circle')
        .data(variants, get_variant_id); // define data-joining 2nd method to allow move-to-front
    window._debug.selection = selection;
    var variant_circles = selection.enter()
        .append('ellipse')
        .attr('class', 'variant-circle')
        .attr('ry', 6)
        .attr('rx', 6)
        .attr('cy', variant_plot_height/2)
        .attr('cx', function(d) { return window.model.plot.x(d.pos); })
        .style('fill', variant_color)
        .style('stroke', 'orange')
        .style('stroke-width', 0)
        .attr('id', get_variant_plot_id)
        .on('mouseover', function(variant) {
            window.model.plot.oval_tip.show(variant);
            variant_plot_default_style();
            variant_plot_hilited_style(d3.select(this));
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
            variant_plot_default_style(d3.select(this));
            window.model.plot.oval_tip.hide(variant);
            window.model.tbl.rows()[0].forEach(function(row_idx) {
                var cur_var = window.model.tbl.row(row_idx).data();
                if (cur_var.pos === variant.pos && cur_var.ref === variant.ref && cur_var.alt === variant.alt) {
                    $(window.model.tbl.row(row_idx).nodes()).removeClass('highlight');
                }
            });
        })
    variant_plot_default_style(variant_circles);
    selection.exit()
        .remove()
    order_selection_by_data_reversed(selection);
}

function order_selection_by_data_reversed(selection){
    // this is like selection.order() but it uses the reverse order.
    // it puts elements from the beginning of the data at the end in the DOM, which is rendered last (on top).
    var n_datasets = selection.length;
    for(var dataset_idx=0;dataset_idx<n_datasets;dataset_idx++) {
        var elems=selection[dataset_idx];
        for(var idx=1;idx<elems.length;idx++) {
            var better = elems[idx-1];
            var lesser = elems[idx];
            if (better && lesser && lesser.nextSibling !== better)
                better.parentNode.insertBefore(lesser, better); // put lesser before better
        }
    }
}

function variant_plot_default_style(selection) {
    if (typeof selection === "undefined") { selection = d3.selectAll('.variant-circle'); }
    return selection
        .style('stroke-width', 0);
}

function variant_plot_hilited_style(selection) {
    return selection
        .style('stroke-width', '3px');
}



function create_variant_table() {
    var columns = [
        {
            title: 'Alleles (rsID)', name: 'allele',
            searchable: false, orderable: false,
            render: (function() {
                var template = _.template(
                    '<a href="<%= window.model.url_prefix %>variant/<%= window.model.chrom %>-<%= variant.pos %>-<%= variant.ref %>-<%= variant.alt %>" target="_blank">'+
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

        // },{
        //     title: 'HGVS', name: 'hgvs',
        //     searchable: false, orderable: false, className: 'dt-head-center',
        //     render: (function() {
        //         var template = _.template(
        //             '<% if (variant.HGVS.length > 20) { %><%= variant.HGVS.slice(0,20) %>...<% } else { %><%= variant.HGVS %><% } %>',
        //             {variable:'variant'});
        //         return function(cell_data, type, row) { return template(row); };
        //     })(),

        },{
            title: 'Consequence', name: 'csq',
            data: 'worst_csqidx', searchable:true, orderable:true, className: 'dt-pad-left',
            render: function(cell_data, type, row) {
                return '<b>' + fmt_annotation(cell_data) + '</b>' +
                    (row.HGVS ? ' ('+row.HGVS+')' : '');
            },

        },{
            title: 'CADD', name:'cadd_phred',
            data: 'cadd_phred', orderable:true, orderSequence:['desc', 'asc'], className: 'dt-right',
            render: function(cell_data, type, row) { return (cell_data===null)?'':cell_data.toFixed(0); },

        },{
            title: 'Quality', name: 'filter',
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

    var update_filter_info = function() {
        window.model.filter_info.maf_ge = parseFloat($('input#maf_ge').val()) / 100; // %
        window.model.filter_info.maf_le = parseFloat($('input#maf_le').val()) / 100; // %
        window.model.filter_info.filter_value = $('select#filter_value').val();
        window.model.filter_info.category = $('#vtf_category > .btn.active').text();
    };
    update_filter_info();

    window.model.tbl = $('#variant_table').DataTable({
        serverSide: true, /* API does all the real work */

        processing: true, /* show "processing" over table while waiting for API */
        deferRender: true, /* only render rows when they're being displayed */

        paging: true,
        pagingType: 'full', /* [first, prev, next, last] */
        pageLength: 100,
        lengthMenu: [10, 100, 1000],

        searching: false,

        scrollX: true,

        ajax: {
            url: window.model.url_prefix + 'api/variants' + window.model.url_suffix,
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

        order: [
            [columns.findIndex(function(d){return d.name==='csq'}), 'asc'],
            [columns.findIndex(function(d){return d.name==='cadd_phred'}), 'desc'],
        ],
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

    $('.variant_table_filter').on('change', function() { update_filter_info(); window.model.tbl.draw(); });

    // hilite corresponding variant-plot circle
    $('#variant_table tbody').on('mouseleave', 'tr', function() {
        var variant = window.model.tbl.row(this).data();
        if (variant) { variant_plot_default_style(d3.select('#' + get_variant_plot_id(variant))); }
        mouse_guide.hide();
    });
    $('#variant_table tbody').on('mouseenter', 'tr', function() {
        var variant = window.model.tbl.row(this).data();
        variant_plot_default_style();
        if (variant) {
            var selection = d3.select('#' + get_variant_plot_id(variant));
            selection.moveToFront();
            variant_plot_hilited_style(selection);
            mouse_guide.show_at(selection.attr('cx'));
        }
    });
}


var mouse_guide = {
    show_at: function(x) {
        d3.selectAll('.mouse_guide').attr('x', x - 1);
        genome_coord = Math.round(window.model.plot.x.invert(x));
        if (x < 0 || x > window.model.plot.genome_coords_width) {
            mouse_guide.hide();
        } else {
            d3.select('#pos_plot_text')
                .attr('transform', fmt('translate({0},{1})', x, pos_plot_height))
                .text(group_thousands(genome_coord));
        }
    },
    hide: function(x) {
        d3.selectAll('.mouse_guide').attr('x', -999);
        d3.select('#pos_plot_text').text('');
    },
    init: function() {
        create_pos_plot();
        d3.selectAll('.genome_g')
            .on('mousemove', function() {var coords = d3.mouse(d3.select(this).node()); mouse_guide.show_at(coords[0]); })
            .on('mouseleave', mouse_guide.hide)
            .insert('rect',':first-child').attr('class', 'genome_g_mouse_catcher'); // recieves mousemove and bubbles it up to `.genome_g`
    },
};




function create_summary_table() {
    var summary_XHR = $.getJSON(window.model.url_prefix + 'api/summary' + window.model.url_suffix);
    $(function() {
        summary_XHR
            .done(function(summary) {
                $('#summary_table').DataTable({
                    paging: false, searching: false, info: false, ordering: false,
                    data: summary,
                    columns: [
                        {title: 'variant type'},
                        {title: 'count (PASS-only)', className:'dt-right', render: function(cell_data, type, row) { return group_thousands_html(cell_data); }}
                    ],
                });
            })
            .fail(function() { console.error('summary XHR failed'); });
    });
}
