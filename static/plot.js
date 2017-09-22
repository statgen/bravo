var genome_coords_margin = {left: 80, right: 5};


function bootstrap_plot() {
    window.model = window.model || {};
    window.model.plot = window.model.plot || {};
    if (!window.model.plot.svg_width) {
        window.model.plot.svg_width = $('#'+coverage_plot.container_id).width();
        window.model.plot.genome_coords_width = window.model.plot.svg_width - genome_coords_margin.left - genome_coords_margin.right;
    }
    if (!window.model.plot.x) {
        // window.model.plot.x() converts from chromosome coordinates -> pixels on the screen.
        var total_range_length = sum(window.model.intervalset.list_of_pairs.map(function(pair){return pair[1]-pair[0]}));
        var SPACING = total_range_length/50; // 2% of the screen will be buffer between intervals
        var domain = [], range = [];
        // first, collapse the space between the intervals
        window.model.intervalset.list_of_pairs.forEach(function(pair) {
            var range_start = (range.length===0 ? 0 : range[range.length-1]+SPACING);
            var range_end = range_start + pair[1]-pair[0];
            domain.push(pair[0]); range.push(range_start);
            domain.push(pair[1]); range.push(range_end);
        });
        // second, convert any values lower or higher than our domain to be just at the edge of the range (to help with exons on region.html, mostly)
        domain = [].concat([-1e100],    domain, [1e100]);
        range =  [].concat([range[0]-1], range, [range[range.length-1]+1]);
        // third, stretch the range to be the width of the screen where we're plotting
        range = range.map(function(x) { return x * window.model.plot.genome_coords_width / range[range.length-1]; });
        window.model.plot.x = d3.scale.linear().domain(domain).range(range);
    }
}


var coverage_plot = {
    container_id: 'coverage_plot_container',
    height: 70,
    margin: {top: 8, bottom: 10},
    color: '#ffa37c',
    create: function() {
        var XHR = $.getJSON(window.model.url_prefix + 'api/coverage' + window.model.url_suffix);
        $(function() {
            bootstrap_plot();

            var svg = d3.select('#'+this.container_id).append("svg")
                .attr("width", window.model.plot.svg_width)
                .attr("height", this.height + this.margin.top + this.margin.bottom)
                .style('display', 'block');
            var genome_g = svg.append("g")
                .attr('id', 'inner_graph')
                .attr('class', 'genome_g')
                .attr("transform", "translate(" + genome_coords_margin.left + "," + this.margin.top + ")");
            genome_g.append('clipPath')
                .attr('id', 'cov-plot-clip')
                .append('rect')
                .attr('x', 0)
                .attr('width', window.model.plot.genome_coords_width)
                .attr('y', 0)
                .attr('height', this.height);

            var data_g = genome_g.append('g')
                .attr('clip-path', 'url(#cov-plot-clip)')
                .attr('id', 'cov-bar-g');

            mouse_guide.register(genome_g);

            var loading_text = genome_g.append('text')
                .attr('text-anchor','middle').text('loading...')
                .attr('transform', fmt('translate({0},{1})', window.model.plot.genome_coords_width/2, this.height/2))

            XHR
                .done(function(coverage_stats) {
                    loading_text.remove();
                    window.model.coverage_stats = coverage_stats;
                    if (window.model.coverage_stats !== null) this.populate(data_g, genome_g);
                    if (window.model.coverage_stats.length === 0) this.no_coverage(genome_g);
                }.bind(this))
                .fail(function() { this.no_coverage(genome_g) });
        }.bind(this));
    },

    no_coverage: function(genome_g) {
        genome_g.append('text')
            .attr('text-anchor', 'middle').text('No Coverage')
            .attr('transform', fmt('translate({0},{1})', window.model.plot.genome_coords_width/2, this.height/2))
    },

    populate: function(data_g, genome_g) {
        var metric = 'mean';
        var max_cov = 1;
        if (metric === 'mean' || metric === 'median') {
            max_cov = d3.max(window.model.coverage_stats, function(d) { return d[metric]; });
        }
        var y = d3.scale.linear()
            .domain([0, max_cov])
            .range([this.height, 0]);
        var yAxis = this.y_axis(y, metric);

        data_g
            .selectAll("rect.cov_plot_bars")
            .data(window.model.coverage_stats)
            .enter()
            .append("rect")
            .attr('class', 'cov_plot_bars')
            .style("fill", this.color)
            .attr("x", function(d) { return window.model.plot.x(d.start); })
            .attr("width", function(d) { return window.model.plot.x(d.end + 1) - window.model.plot.x(d.start) + 1; })
            .attr("y", function(d) { return y(d[metric]) || 0; })
            .attr("height", function(d) { return (this.height - y(d[metric])) || 0; }.bind(this));

        genome_g.append("g")
            .attr("class", "y axis")
            .call(yAxis);

        // Handle changes
        $('.coverage_metric_buttons').change(function () {
            coverage_plot.change_metric_type($(this).attr('id').replace('_covmet_button', ''));
        });
        $('#over_x_select').change(function () {
            coverage_plot.change_metric($(this).val().replace('X', ''));
        });
        $('#average_select').change(function () {
            coverage_plot.change_metric($(this).val());
        });
    },
    y_axis: function(y_scale, metric) {
        var yAxis = d3.svg.axis()
            .scale(y_scale)
            .orient('left')
            .ticks(3);
        if (metric === 'mean' || metric === 'median')
            yAxis = yAxis.tickFormat(function(d) {return d.toString() + '\u00d7'});
        else
            yAxis = yAxis.tickFormat(d3.format('%'));
        return yAxis;
    },

    change_metric_type: function(metric_type) {
        $('.coverage_subcat_selectors').css('display', 'none');
        if (metric_type == 'covered') {
            $('#over_x_select_container').css('display', 'inline-block');
            var metric = $('#over_x_select').val().replace('X', '');
        } else {
            $('#average_select_container').css('display', 'inline-block');
            var metric = $("#average_select").val();
        }
        this.change_metric(metric);
    },
    change_metric: function(metric) {
        var max_cov = 1;
        if (metric === 'mean' || metric === 'median') {
            max_cov = d3.max(window.model.coverage_stats, function(d) { return d[metric]; });
        }

        var y = d3.scale.linear()
            .domain([0, max_cov])
            .range([this.height, 0]);

        var svg = d3.select('#'+this.container_id).select('#inner_graph');

        var yAxis = this.y_axis(y, metric);
        svg.select(".y.axis").call(yAxis);

        svg.select('#cov-bar-g')
            .selectAll("rect.cov_plot_bars")
            .data(window.model.coverage_stats)
            .attr("y", function(d) { return y(d[metric]); })
            .attr("height", function(d) { return this.height - y(d[metric]); }.bind(this));
    },
};


var transcripts_plot = {
    height: 15, // not including margins
    margin: {top:1, bottom:1},
    label_gradient_width: 40,
    colors: {canonical:'darkblue', coding:'steelblue', noncoding:'lightsteelblue'},
    create: function() {
        bootstrap_plot();

        var num_transcripts = this.count_transcripts(window.model.genes);
        var gotta_use_subset = num_transcripts > 7;
        var genes = (gotta_use_subset ? this.get_genes_subset() : window.model.genes)
        this.render_genes(genes);
        if (gotta_use_subset) {
            this.add_draw_all_button(num_transcripts - this.count_transcripts(genes), window.model.genes.length - genes.length);
        }
    },
    count_transcripts: function(genes){
        return sum(genes.map(function(gene) {return gene.transcripts.length}));
    },
    get_genes_subset: function() {
        // if there are fewer than 5 transcripts we'll probably have an infinite loop so don't do that.
        var genes = window.model.genes.slice(0,5).map(function(gene) {
            gene = deepcopy(gene);
            gene.transcripts = [gene.transcripts[0]];
            return gene;
        });
        for (var trans_i=1; ; trans_i++) {
            for (var gene_i=0; gene_i < window.model.genes.length; gene_i++) {
                if (this.count_transcripts(genes) >= 5) return genes;
                var trans = window.model.genes[gene_i].transcripts[trans_i];
                if (trans) genes[gene_i].transcripts.push(trans);
            }
        }
    },
    render_genes: function(genes) {
        genes.forEach(function(gene) {
            gene.transcripts.forEach(function(transcript, i) {
                var is_coding = _.any(transcript.exons, function(exon) {return exon.feature_type === 'CDS'});
                this.create_one(gene, transcript, is_coding);
            }.bind(this))
        }.bind(this))
    },
    create_one: function(gene, transcript, is_coding) {
        var color = (transcript.canonical ? this.colors.canonical : (is_coding ? this.colors.coding : this.colors.noncoding));

        var div = d3.select('#transcripts_plot_container').append('div');
        var label_p = div.append('p').attr('class','transcript-label').style('position','absolute').style('margin',0);
        var svg = div.append('svg')
            .attr('width', window.model.plot.svg_width)
            .attr('height', this.height + this.margin.top + this.margin.bottom)
            .style('display', 'block')
            .style('position','relative').style('z-index',1) // allows z-index
            .style('pointer-events','none') // don't capture mouseover, to let it pass thru to label
        svg.append('rect')
            .attr('transform', fmt('translate({0},0)', genome_coords_margin.left))
            .style('fill','white').attr('height','100%').attr('width','100%')
            .style('pointer-events','all');
        this._populate_label(gene, transcript, is_coding, label_p, svg, color);

        var genome_g = svg.append('g')
            .attr('id', 'gene_track')
            .attr('class', 'genome_g')
            .attr('transform', fmt('translate({0},{1})', genome_coords_margin.left, this.margin.top))
            .style('pointer-events', 'all');
        var data_g = genome_g.append('g');
        mouse_guide.register(genome_g);

        var exon_tip = d3.tip()
            .attr('class', 'd3-tip')
            .style('z-index', 2)
            .html(function(d) {
                return transcript.transcript_id + '<br>' +
                    (d.feature_type==='CDS'?'Coding Sequence':d.feature_type) + '<br>' +
                    'start: ' + group_thousands_html(d.start) + '<br>' +
                    'stop: ' + group_thousands_html(d.stop) + '<br>' +
                    'strand: ' + d.strand;
            });
        svg.call(exon_tip);

        var intervals_to_show = window.model.intervalset.list_of_pairs
            .filter(function(pair) { return pair[0] <= transcript.stop && pair[1] >= transcript.start })
            .map(function(pair) { return [Math.max(pair[0], transcript.start), Math.min(pair[1], transcript.stop)] });
        data_g.selectAll('line.intervals')
            .data(intervals_to_show)
            .enter()
            .append('line')
            .attr("y1", this.height/2)
            .attr("y2", this.height/2)
            .attr("x1", function(d) { return window.model.plot.x(d[0]) })
            .attr("x2", function(d) { return window.model.plot.x(d[1]) })
            .attr("stroke-width", 1)
            .attr("stroke", color);

        data_g.selectAll('rect.exon')
            .data(transcript.exons)
            .enter()
            .append('a')
            .attr('xlink:href', function(d){return fmt('{0}region/{1}-{2}-{3}', window.model.url_prefix, window.model.intervalset.chrom, d.start, d.stop)})
            .append('rect')
            .attr('class', 'exon')
            .style('fill', color)
            .attr('y', function(d){return d.feature_type==='CDS' ? 0 : this.height/4}.bind(this))
            .attr('height',function(d){return d.feature_type==='CDS' ? this.height : this.height/2}.bind(this))
            .each(function(d) {
                var start = window.model.plot.x(d.start);
                var width = window.model.plot.x(d.stop) - start;
                d3.select(this).attr('x', start).attr('width', width)
            })
            .on('mouseover', exon_tip.show)
            .on('mouseout', exon_tip.hide)
    },
    _populate_label: function(gene, transcript, is_coding, label_p, svg, color) {
        var grad = svg.append('linearGradient').attr('id', 'label-mask-gradient').attr('x2', 1)
        grad.append('stop').attr('offset',0).attr('stop-opacity',0).attr('stop-color','white')
        grad.append('stop').attr('offset',0.8).attr('stop-opacity',1).attr('stop-color','white')
        grad.append('stop').attr('offset',1).attr('stop-opacity',1).attr('stop-color','white')
        svg.append('rect')
            .attr('x', genome_coords_margin.left-this.label_gradient_width)
            .attr('width', this.label_gradient_width)
            .attr('height', '100%')
            .attr('fill', 'url(#label-mask-gradient)')

        var font_weight = transcript.canonical ? 'bold' : 'inherit';
        label_p
            .style('height', fmt('{0}px', this.height + this.margin.bottom + this.margin.top))
            .on('mouseover', function() { d3.selectAll('.transcript-label').style('z-index', 2) })
            .on('mouseout', function() { d3.selectAll('.transcript-label').style('z-index', 0) })
            .style('background-color','white')
            .style('padding-right', '0.5em')
        var show_gene = window.model.url_suffix.startsWith('/region/'); // TODO: use a cleaner method to check whether we're on region.html
        if (show_gene) {
            label_p.append('a')
                .attr('href', fmt('{0}gene/{1}', window.model.url_prefix, gene.gene_id))
                .style('font-style','italic')
                .style('margin-right', '0.2em')
                .text(gene.gene_name || gene.gene_id)
                .style('color', color)
                .style('font-weight', font_weight)
        }
        var label = label_p.append('a')
            .attr('href', fmt('{0}transcript/{1}', window.model.url_prefix, transcript.transcript_id))
            .text(transcript.transcript_id + (transcript.canonical?' (canonical)':''))
            .style('color', color)
            .style('font-weight', font_weight)
    },
    add_draw_all_button: function(num_transcripts_remaining, num_genes_remaining) {
        var text = fmt('Show remaining {0} transcripts',num_transcripts_remaining) +
            (num_genes_remaining ? fmt(' on {0} more gene{1}',num_genes_remaining,(num_genes_remaining===1?'':'s')) : '');

        d3.select('#transcripts_plot_container')
            .append('div')
            .style('text-align', 'center')
            .append('button')
            .attr('class', 'btn btn-default')
            .text(text)
            .on('click', function(){
                d3.select('#transcripts_plot_container').selectAll('*').remove()
                this.render_genes(window.model.genes);
            }.bind(this))
    },
};


var pos_plot = {
    container_id: 'pos_plot_container',
    height: 15,
    margin: {top: 0, bottom: 0},
    create: function() {
        bootstrap_plot();
        var svg = d3.select('#'+this.container_id).append('svg')
            .attr('width', window.model.plot.svg_width)
            .attr('height', this.height + this.margin.top + this.margin.bottom)
            .style('display', 'block');
        var genome_g = svg.append('g')
            .attr('class', 'genome_g')
            .attr('transform', fmt('translate({0},{1})', genome_coords_margin.left, this.margin.top));
        genome_g.append('text').attr('id', 'pos_plot_text')
            .attr('text-anchor', 'middle');
        mouse_guide.register(genome_g, true);
    }
};



var variant_plot = {
    container_id: 'variant_plot_container',
    height: 20,
    margin: {top: 0, bottom: 0},

    create: function() {
        bootstrap_plot();

        var svg = d3.select('#'+this.container_id).append("svg")
            .attr("width", window.model.plot.svg_width)
            .attr("height", this.height + this.margin.top + this.margin.bottom)
            .style('display', 'block');
        var genome_g = svg.append("g")
            .attr('id', 'variant_track')
            .attr('class', 'genome_g')
            .attr("transform", fmt('translate({0},{1})', genome_coords_margin.left, this.margin.top));

        genome_g.append('clipPath')
            .attr('id', 'variant-plot-clip')
            .append('rect')
            .attr('x', 0)
            .attr('width', window.model.plot.genome_coords_width)
            .attr('y', 0)
            .attr('height', this.height);

        mouse_guide.register(genome_g);

        genome_g.selectAll('line.intervals')
            .data(window.model.intervalset.list_of_pairs)
            .enter()
            .append('line')
            .attr("y1", this.height/2)
            .attr("y2", this.height/2)
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
        });
        svg.call(window.model.plot.oval_tip);
    },
    get_variant_name: function(variant) { return ''+variant.pos+'-'+variant.ref+'-'+variant.alt },
    get_variant_id: function(variant) { return 'variant-plot-'+this.get_variant_name(variant) },
    change: function(variants) {
        var selection = d3.select('#variant_track .data_g')
            .selectAll('.variant-circle')
            .data(variants, function(variant) { return this.get_variant_name(variant) }.bind(this)); // define data-joining 2nd method to allow move-to-front
        var variant_circles = selection.enter()
            .append('ellipse')
            .attr('class', 'variant-circle')
            .attr('ry', 6)
            .attr('rx', 6)
            .attr('cy', this.height/2)
            .attr('cx', function(d) { return window.model.plot.x(d.pos); })
            .style('fill', variant_color)
            .style('stroke', 'orange')
            .style('stroke-width', 0)
            .attr('id', function(d) { return this.get_variant_id(d) }.bind(this))
            .on('mouseover', function(variant) {
                window.model.plot.oval_tip.show(variant);
                variant_plot.default_style();
                variant_plot.hilited_style(d3.select(this));
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
                variant_plot.default_style(d3.select(this));
                window.model.plot.oval_tip.hide(variant);
                window.model.tbl.rows()[0].forEach(function(row_idx) {
                    var cur_var = window.model.tbl.row(row_idx).data();
                    if (cur_var.pos === variant.pos && cur_var.ref === variant.ref && cur_var.alt === variant.alt) {
                        $(window.model.tbl.row(row_idx).nodes()).removeClass('highlight');
                    }
                });
            })
        this.default_style(variant_circles);
        selection.exit().remove()
        this.order_selection_by_data_reversed(selection);
    },

    default_style: function(selection) {
        return (selection || d3.selectAll('.variant-circle')).style('stroke-width', 0);
    },
    hilited_style: function(selection) {
        return selection.style('stroke-width', '3px');
    },

    order_selection_by_data_reversed: function(selection){
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
    },
};


var variant_table = {
    columns: [
        {
            title: 'Alleles (rsID)', name: 'allele',
            searchable: false, orderable: false,
            render: (function() {
                var template = _.template(
                    '<a href="<%= window.model.url_prefix %>variant/<%= window.model.intervalset.chrom %>-<%= variant.pos %>-<%= variant.ref %>-<%= variant.alt %>" target="_blank">'+
                        '<% if (variant.ref.length > 20) { %><%= variant.ref.slice(0,20) %>...<% } else { %><%= variant.ref %><% } %> / '+
                        '<% if (variant.alt.length > 20) { %><%= variant.alt.slice(0,20) %>...<% } else { %><%= variant.alt %><% } %>'+
                        '<% if (variant.rsids.length) { %> (<%= variant.rsids.join(", ") %>)<% } %>'+
                        '</a>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Position (chr' + window.model.intervalset.chrom + ')', name: 'pos',
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
                    (row.HGVS ? ' ('+row.HGVS+')' : '') +
                    (row.low_conf ? ' <i>(Low&nbsp;Confidence)</i>':'');
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
            title: 'Frequency (%)', name: 'allele_freq',
            data: 'allele_freq', searchable:true, orderable:true, orderSequence:['desc','asc'], className: 'dt-pad-left',
            render: function(cell_data, type, row) { return perc_sigfigs_html(cell_data, 2); },

        },
    ],

    update_filter_info: function() {
        window.model.filter_info.maf_ge = parseFloat($('input#maf_ge').val()) / 100; // %
        window.model.filter_info.maf_le = parseFloat($('input#maf_le').val()) / 100; // %
        window.model.filter_info.filter_value = $('select#filter_value').val();
        window.model.filter_info.category = $('#vtf_category > .btn.active').text();
    },
    create: function() {

        window.model = window.model || {};
        window.model.filter_info = window.model.filter_info || {};

        this.update_filter_info();

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
                    variant_plot.change(resp.data);
                    return resp.data;
                },
            },

            order: [
                [this.columns.findIndex(function(d){return d.name==='csq'}), 'asc'],
                [this.columns.findIndex(function(d){return d.name==='cadd_phred'}), 'desc'],
            ],
            columns: this.columns,

            dom: '<ipl>rft', // default is 'lfrtip'.  l=length f=filtering t=table i=info p=paging, r=processing

            language: {
                info: 'Showing variants _START_ - _END_ of _TOTAL_',
                infoFiltered: '(filtered from _MAX_ variants)',
                infoEmpty: 'No matching variants',
                thousands: '\u202f',
                lengthMenu: 'Show _MENU_ variants',
            }

        });

        $('.variant_table_filter').on('change', function() { variant_table.update_filter_info(); window.model.tbl.draw(); });

        // hilite corresponding variant-plot circle
        $('#variant_table tbody').on('mouseleave', 'tr', function() {
            var variant = window.model.tbl.row(this).data();
            if (variant) { variant_plot.default_style(d3.select('#' + variant_plot.get_variant_id(variant))); }
            mouse_guide.hide();
        });
        $('#variant_table tbody').on('mouseenter', 'tr', function() {
            var variant = window.model.tbl.row(this).data();
            variant_plot.default_style();
            if (variant) {
                var selection = d3.select('#' + variant_plot.get_variant_id(variant));
                selection.moveToFront();
                variant_plot.hilited_style(selection);
                mouse_guide.show_at(selection.attr('cx'));
            }
        });
    },
};


var mouse_guide = {
    register: function(genome_g, no_line) {
        if (!no_line) {
            genome_g
                .append('rect').attr('class', 'mouse_guide')
                .attr('x', -999).attr('width','2px').attr('height','100%')
        }
        genome_g
            .on('mousemove', function() {var coords = d3.mouse(d3.select(this).node()); mouse_guide.show_at(coords[0]); })
            .on('mouseleave', mouse_guide.hide)
            .insert('rect',':first-child').attr('class', 'genome_g_mouse_catcher'); // recieves mousemove and bubbles it up to `.genome_g`
    },
    show_at: function(x) {
        if (this._x_is_in_intervalset(x)) {
            d3.selectAll('.mouse_guide').attr('x', x - 1);
            var genome_coord = Math.round(window.model.plot.x.invert(x));
            d3.select('#pos_plot_text')
                .attr('transform', fmt('translate({0},{1})', x, pos_plot.height))
                .text(group_thousands(genome_coord));
        } else {
            mouse_guide.hide();
        }
    },
    hide: function(x) {
        d3.selectAll('.mouse_guide').attr('x', -999);
        d3.select('#pos_plot_text').text('');
    },
    _x_is_in_intervalset: function(x) {
        // TODO: find a way to make this function faster.  Maybe bisection?
        if (x < 0 || x > window.model.plot.genome_coords_width) { return false; }
        var pairs = window.model.intervalset.list_of_pairs;
        var genome_coord = Math.round(window.model.plot.x.invert(x));
        return pairs[0][0] <= genome_coord &&
            genome_coord <= pairs[pairs.length-1][1] &&
            _.any(pairs, pair => (pair[0] <= genome_coord && genome_coord <= pair[1]));
    }
}


var summary_table = {
    container_id: 'summary_table',
    populate: function(summary) {
        $('#'+summary_table.container_id).DataTable({
            paging: false, searching: false, info: false, ordering: false,
            data: summary,
            columns: [
                {title: 'variant type'},
                {title: 'count (PASS-only)', className:'dt-right', render: function(cell_data, type, row) { return group_thousands_html(cell_data); }}
            ],
        });
    },
};


var density_plot = {
    container_id: 'density_plot_container',
    height: 60,
    margin: {top: 5, bottom: 15},
    color: {variants:'#999', singletons:'#ddd', '.1%+':'#555'},
    create: function() {
        if (window.model.intervalset.list_of_pairs.length === 1 && window.model.intervalset.list_of_pairs[0][1] - window.model.intervalset.list_of_pairs[0][0] < 200) return;
        bootstrap_plot();
        var svg = d3.select('#'+density_plot.container_id).append('svg')
            .attr("width", window.model.plot.svg_width)
            .attr("height", this.height + this.margin.top + this.margin.bottom)
            .style('display', 'block');
        var genome_g = svg.append("g")
            .attr('class', 'genome_g')
            .attr("transform", "translate(" + genome_coords_margin.left + "," + this.margin.top + ")");
        var data_g = genome_g.append('g')
            .attr('class', 'data_g');
        mouse_guide.register(genome_g);

        genome_g.append('text')
            .attr('class', 'loading')
            .attr('text-anchor','middle').text('loading...')
            .attr('transform', fmt('translate({0},{1})', window.model.plot.genome_coords_width/2, this.height/2))

    },
    populate: function(interval_summaries) {
        interval_summaries.forEach(function(isumm) {
            isumm.length = isumm.stop - isumm.start + 1;
            isumm.variant_density = isumm.variants / isumm.length * 1000;
            isumm.singleton_density = isumm.singletons / isumm.length * 1000;
            isumm['.1%+_density'] = isumm['.1%+'] / isumm.length * 1000;
        });
        var svg = d3.select('#'+density_plot.container_id).select('svg');
        if (svg.empty()) return;
        var genome_g = svg.select('.genome_g');
        var data_g = svg.select('.data_g');
        genome_g.select('.loading').remove();
        var tip = d3.tip()
            .attr('class', 'd3-tip')
            .html(function(isumm) {
                return 'start: ' + group_thousands_html(isumm.start) + '<br>' +
                    'stop: ' + group_thousands_html(isumm.stop) + '<br>' +
                    group_thousands_html(isumm.length) + ' bases<br><br>' +
                    group_thousands_html(isumm.variants) + ' variants<br>' +
                    group_thousands_html(isumm.singletons) + ' singletons<br>' +
                    group_thousands_html(isumm['.1%+']) + ' AF>0.1%<br><br>' +
                    group_thousands_html(Math.round(isumm.variant_density)) + ' variants/kb<br>' +
                    group_thousands_html(Math.round(isumm.singleton_density)) + ' singletons/kb<br>' +
                    group_thousands_html(Math.round(isumm['.1%+_density'])) + ' AF>0.1%/kb<br><br>' +
                    (isumm.variants ? perc_sigfigs_html(isumm.singletons/isumm.variants, 2, 0) + ' singletons' : '')
            });
        svg.call(tip);

        var y = d3.scale.linear()
            .domain([0, d3.max(interval_summaries, function(isumm){return isumm.variant_density})])
            .range([this.height, 0]);
        data_g
            .selectAll('rect.tip_target')
            .data(interval_summaries)
            .enter()
            .append('rect')
            .attr('class', 'tip_target')
            .style('fill', 'white')
            .attr('x', function(d) { return window.model.plot.x(d.start) })
            .attr('width', function(d) { return window.model.plot.x(d.stop + 1) - window.model.plot.x(d.start) + 1; })
            .attr('y', function(d) { return 0 })
            .attr('height', function(d) { return density_plot.height })
            .on('mouseover', tip.show)
            .on('mouseout', tip.hide);
        data_g
            .selectAll('rect.num_variants')
            .data(interval_summaries)
            .enter()
            .append('rect')
            .attr('class', 'num_variants')
            .style('fill', density_plot.color.variants)
            .attr('x', function(d) { return window.model.plot.x(d.start) })
            .attr('width', function(d) { return window.model.plot.x(d.stop + 1) - window.model.plot.x(d.start) + 1; })
            .attr('y', function(d) { return y(d.variant_density) })
            .attr('height', function(d) { return (density_plot.height - y(d.variant_density)) || 0 })
            .style('pointer-events', 'none');
        data_g
            .selectAll('rect.num_singletons')
            .data(interval_summaries)
            .enter()
            .append('rect')
            .attr('class', 'num_singletons')
            .style('fill', density_plot.color.singletons)
            .attr('x', function(d) { return window.model.plot.x(d.start) })
            .attr('width', function(d) { return window.model.plot.x(d.stop + 1) - window.model.plot.x(d.start) + 1; })
            .attr('y', function(d) { return y(d.singleton_density) })
            .attr('height', function(d) { return (density_plot.height - y(d.singleton_density)) || 0 })
            .style('pointer-events', 'none');
        data_g
            .selectAll('rect.num_common')
            .data(interval_summaries)
            .enter()
            .append('rect')
            .attr('class', 'num_common')
            .style('fill', density_plot.color['.1%+'])
            .attr('x', function(d) { return window.model.plot.x(d.start) })
            .attr('width', function(d) { return window.model.plot.x(d.stop + 1) - window.model.plot.x(d.start) + 1; })
            .attr('y', function(d) { return y(d.variant_density) })
            .attr('height', function(d) { return (density_plot.height - y(d['.1%+_density'])) || 0 })
            .style('pointer-events', 'none');

        var yAxis = d3.svg.axis().scale(y).orient('left').ticks(2).tickFormat(function(d){return d+'/kb'});
        genome_g.append('g')
            .attr('class', 'y axis')
            .call(yAxis)
    }
};


var summary_XHR = $.getJSON(window.model.url_prefix + 'api/summary' + window.model.url_suffix);
$(function() {
    density_plot.create();
    summary_XHR
        .done(function(summary_resp) {
            summary_table.populate(summary_resp.summary);
            density_plot.populate(summary_resp.interval_summaries);
            // singleton_density_plot.populate(summary_resp.interval_summaries);
        })
        .fail(function() { console.error('summary XHR failed'); });
});
coverage_plot.create();
$(function() {
    transcripts_plot.create();
    variant_plot.create();
    variant_table.create();
    pos_plot.create();
});
