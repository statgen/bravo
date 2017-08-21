var gene_chart_margin = {top: 10, right: 30, bottom: 30, left: 80};
var gene_chart_margin_lower = {top: 5, right: gene_chart_margin.right, bottom: 5, left: gene_chart_margin.left};
var lower_gene_chart_height = 50 - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom;
var gene_chart_height = 300 - gene_chart_margin.top - gene_chart_margin.bottom - lower_gene_chart_height - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom;

function bootstrap_plot() {
    window.model = window.model || {};
    window.model.plot = window.model.plot || {};
    if (!window.model.plot.width) {
        var width_including_margin = $('#region_plot_container').width();
        window.model.plot.width = width_including_margin - gene_chart_margin.left - gene_chart_margin.right;
    }
    if (!window.model.plot.x) {
        window.model.plot.x = d3.scale.linear()
            .domain([window.model.start, window.model.stop])
            .range([0, window.model.plot.width]);
    }
}

function create_region_chart(cov_data) {
    bootstrap_plot();

    var metric = 'mean';

    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(cov_data, function(d) { return d[metric]; });
    }

    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([gene_chart_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(window.model.plot.x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#region_plot_container').append("svg")
        .attr("width", window.model.plot.width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + gene_chart_margin.left + "," + gene_chart_margin.top + ")");

    svg.selectAll("bar")
        .data(cov_data)
        .enter()
        .append("rect")
        .attr('class', 'main_plot_bars')
        .style("fill", "steelblue")
        .attr("x", function(d) {
            return window.model.plot.x(d.start);
        })
        .attr("width", function(d) {
            var length_in_bases = d.end - d.start + 1;
            var width_of_base = window.model.plot.width/cov_data.length;
            return length_in_bases * width_of_base;
        })
        .attr("y", function(d) {
            return y(d[metric]) || 0;
        })
        .attr("height", function(d) {
            return (gene_chart_height - y(d[metric])) || 0;
        });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + gene_chart_height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);
}

function change_region_chart_metric(cov_data, metric) {

    console.log(['cov_data', cov_data]);

    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(cov_data, function(d) { return d[metric]; });
    }
    console.log(max_cov);

    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([gene_chart_height, 0]);

    var svg = d3.select('#region_plot_container').select('#inner_graph');

    svg.selectAll("rect")
        .data(cov_data)
        .transition()
        .duration(500)
        .attr("y", function(d) { return y(d[metric]); })
        .attr("height", function(d) { return gene_chart_height - y(d[metric]); });

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    svg.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);

}

function create_variant_oval_track() {
    bootstrap_plot();

    var svg_outer = d3.select('#region_plot_container').append("svg")
        .attr("width", window.model.plot.width + gene_chart_margin_lower.left + gene_chart_margin_lower.right)
        .attr("height", lower_gene_chart_height)
        .append("g")
        .attr('id', 'track')
        .attr("transform", "translate(" + gene_chart_margin_lower.left + "," + 0 + ")");

    svg_outer.append("line")
        .attr("y1", lower_gene_chart_height/2)
        .attr("y2", lower_gene_chart_height/2)
        .attr("x1", 0)
        .attr("x2", window.model.plot.width)
        .attr("stroke-width", 1)
        .attr("stroke", "lightsteelblue");

    window.model.plot.oval_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        var csq = d.major_consequence.replace(/_/g, ' ');
        if (csq.length > 15) { csq = csq.substr(0, 15) + '...'; } // because d3-tip tooltips fall off the page
        return group_thousands(d.pos) + '<br>' +
            csq + '<br>' +
            (d.filter === 'PASS' ? '' : d.filter + '<br>') +
            'MAF: ' + perc_sigfigs(d.allele_freq, 2);
        //return JSON.stringify(d);
    });
    svg_outer.call(window.model.plot.oval_tip);
}

function change_variant_oval_track(variants) {
    console.log('plotting', variants);
    var selection = d3.select('#track').selectAll('.variant-circle').data(variants);
    selection.enter() //ENTER
        .append('ellipse')
        //.attr('foo', function(d) { console.log('enter-each', d)})
        .attr('class', 'variant-circle')
        .style('opacity', 0.5)
        .style('fill', 'blue')
        .attr('ry', 4)
        .attr('rx', 4)
        .attr('cy', lower_gene_chart_height/2)
        .attr('cx', function(d) { return window.model.plot.x(d.pos); })
        .on('mouseover', window.model.plot.oval_tip.show)
        .on('mouseout', window.model.plot.oval_tip.hide)
    selection // UPDATE // todo: learn how to get union of ENTER+UPDATE
        //.attr('foo', function(d) { console.log('update-each', d)})
        .attr('cx', function(d) { return window.model.plot.x(d.pos); })
    selection.exit() //EXIT
        //.attr('foo', function(d) { console.log('exit-each', d)})
        .remove()
}


$(document).ready(function() {
    if (window.model.coverage_stats != null) {
        create_region_chart(window.model.coverage_stats);
        // Change coverage plot
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
            change_region_chart_metric(window.model.coverage_stats, v);
        });
        $('#over_x_select').change(function () {
            change_region_chart_metric(window.model.coverage_stats, $(this).val().replace('X', ''));
        });
        $('#average_select').change(function () {
            change_region_chart_metric(window.model.coverage_stats, $(this).val());
        });
    }

    create_variant_oval_track();
});




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
                        '<% if (variant.rsids.length) { %>(<%= variant.rsids.join(", ") %>)<% } %>'+
                        '</a>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Position on chr' + window.model.chrom, name: 'pos',
            data: 'pos', searchable: true,  orderable: true, className: 'dt-right',
            render: function(cell_data, type, row) { return group_thousands(cell_data); },

        },{
            title: 'HGVS', name: 'hgvs',
            searchable: false, orderable: false, className: 'dt-head-center',
            render: (function() {
                var template = _.template(
                    '<% if (variant.HGVS != "") { %>' +
                    '<% if (variant.HGVS.length > 20) { %><%= variant.HGVS.slice(0,20) %>...<% } else { %><%= variant.HGVS %><% } %>' +
                        '<% if (variant.CANONICAL !== "YES") { %><span class="tooltip-table-header" data-tooltip="Annotation is for non-canonical transcript">&dagger;</span><% } %>'+
                        '<% } %>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Annotation', name: 'csq',
            data: 'major_consequence', searchable:true, orderable:false, className: 'dt-center',
            render: function(cell_data, type, row) {
                return fmt_annotation(cell_data);
            },

        },{
            title: 'CADD', name:'cadd',
            orderable:false, className: 'dt-right',
            render: function() { return '0'; },

        },{
            title: 'QC', name: 'filter',
            data: 'filter', searchable:true, orderable:false, className: 'dt-center',

        },{
            title: 'HomRef', name: 'n_homref',
            searchable:true, orderable:false, className: 'dt-right',
            render: function(cell_data, type, row) {
                var num_het_samples = row.allele_count - 2*row.hom_count;
                var num_homref_samples = ((row.allele_num - row.allele_count) - num_het_samples) / 2;
                return group_thousands(num_homref_samples);
            },

        },{
            title: 'Het', name: 'n_het',
            searchable:true, orderable:false, className: 'dt-right',
            render: function(cell_data, type, row) {
                var num_het_samples = row.allele_count - 2*row.hom_count;
                return group_thousands(num_het_samples);
            },

        },{
            title: 'HomAlt', name: 'hom_count',
            searchable:true, orderable:true, className: 'dt-right',
            render: function(cell_data, type, row) {
                return group_thousands(row.hom_count);
            },

        },{
            title: 'Frequency', name: 'allele_freq',
            data: 'allele_freq', searchable:true, orderable:true, className: 'dt-pad-left',
            render: function(cell_data, type, row) { return perc_sigfigs(cell_data, 2); },

        },
    ];

    window.model = window.model || {};
    window.model.filter_info = window.model.filter_info || {};
    window.model.filter_info.start = window.model.start;
    window.model.filter_info.stop = window.model.stop;
    window.model.filter_info.chrom = window.model.chrom;

    $('#variant_table').DataTable({
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
                change_variant_oval_track(resp.data);
                return resp.data;
            },
        },

        order: [[columns.map(function(d){return d.data==='pos'}).indexOf(true), 'asc']],
        columns: columns,

        dom: '<ip>ftr'  //'ipftr',  // default is 'lfrtip'.  l=length f=filtering t=table i=info p=paging, r=processing

    });

    $('.variant_table_filter').on('change', function() {
        window.model.filter_info.pos_ge = parseInt($('input#pos_ge').val());
        window.model.filter_info.pos_le = parseInt($('input#pos_le').val());
        window.model.filter_info.filter_value = $('select#filter_value').val();
        console.log(window.model);
        $('#variant_table').DataTable().draw();
    });
}
$(create_variant_table);



function group_thousands(x) {
    // group_thousands_space(1234567.7654321) -> "1 234 567.765 432 1"
    try {
        var thinspace = '\u2009';
        var parts = (''+x).split('.');
        var b = parts[0]; // before the decimal
        var neg = b[0] === '-';
        if (neg) b = b.substr(1);
        var L = b.length;
        var i = b.length;
        var bf = ''; // before the decimal, formatted
        while (i--) bf = (i===0?'':((L-i)%3?'':thinspace)) + b.charAt(i) + bf;
        if (neg) bf = '-' + bf;
        if (parts.length === 1) return bf;
        var a = parts[1]; // after the decimal
        L = a.length;
        i = 0;
        var af = ''; // after the decimal, formatted
        while (i < L) { af += (i===0?'':(i%3?'':thinspace)) + a.charAt(i); i++; }
        return bf + '.' + af;
    } catch(e) {
        console.error(['exception occurrend in group_thousands', e, x]);
        return x;
    }
}
function perc_sigfigs(d, n_sigfigs, n_left_of_decimal) {
    /* like proportion_sigfigs but outputs a percentage.
       left-pads with &nbsp; to line up `100%` with `0%`
       a leading 9 will never be rounded up, so the number of digits can never change. */
    try {
        if (typeof n_sigfigs === 'undefined') n_sigfigs = 2;
        if (typeof n_left_of_decimal === 'undefined') n_left_of_decimal = 3;
        var nobreakspace = '\u2007'; // '\u00A0';

        d = d*100;
        var parts = d.toString().replace(/^0+/, '').split('.');
        var decimal_offset = parts[0].length;
        if (parts[0].length > n_sigfigs) n_sigfigs = parts[0].length; /* I hate zero used to imply "dunno" */
        var onestring = parts.join('');
        onestring = sigfig(onestring, n_sigfigs);
        if (onestring.length < decimal_offset) { throw 'wat ' + onestring + ' ' + decimal_offset }
        var before_decimal = onestring.substr(0, decimal_offset);
        var after_decimal = onestring.substr(decimal_offset);
        if (before_decimal.length === 0) before_decimal = '0';
        before_decimal = leftpad(before_decimal, n_left_of_decimal, nobreakspace);
        var full = before_decimal + (after_decimal.length ? '.' + after_decimal : '');
        full = group_thousands(full);
        return full + '%'; // + ' (' + d.toString() + ')';
    } catch(e) {
        console.error(['exception occurrend in perc_sigfigs', e, [d, n_sigfigs, n_left_of_decimal]]);
        return d;
    }
}
function sigfig(numstring, k) {
    /* sigfigs("999499", 1) -> "9995" */
    var not_sig = ''
    if (numstring[0] === '9') not_sig = numstring.match(/^9*/)[0];
    if (numstring[0] === '0') not_sig = numstring.match(/^0*/)[0];
    var sig = numstring.substr(not_sig.length);
    var firstk = sig.substr(0, k);
    if (sig.length > k) { /* round */
        var rounder = sig[k];
        if (parseInt(rounder) >= 5) firstk = (parseInt(firstk) + 1).toString();
    }
    return not_sig + firstk;
}
function leftpad(str, len, padding) {
    str = String(str); len = len>>0; if (str.length >= len) return str; return padding.repeat(len - str.length) + str;
}
function rightpad(str, len, padding) {
    str = String(str); len = len>>0; if (str.length >= len) return str; return str + padding.repeat(len - str.length);
}
function fmt_annotation(anno) {
    return anno.replace('_variant', '').replace(/_/g, ' ').replace('utr', 'UTR').replace('3 prime', "3'").replace('5 prime', "5'").replace('nc ', "non-coding ");
}
