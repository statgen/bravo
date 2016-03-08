
quality_chart_margin = {top: 10, right: 30, bottom: 45, left: 65};
quality_chart_height = 250 - quality_chart_margin.top - quality_chart_margin.bottom;
quality_chart_width = 300 - quality_chart_margin.left - quality_chart_margin.right;
xoffset = 40;
yoffset = 55;

function draw_region_coverage(raw_data, metric, ref) {
    // pjvh removed a bunch of functionality from this.  There's some useless code left behind.
    // If this function gets a single base, it draws the full distribution.
    // If it receives multiple bases, it draws a coverage graph, letting the user select mean, median, % > X
    // TODO: always draw coverage for only the first base.
    region_chart_width = 500;
    region_chart_margin = {top: 10, right: 50, bottom: 55, left: 50};
    if (raw_data.length > 1) {
        var data = raw_data;
        var chart_width = _.min([region_chart_width, data.length*30]);
        var x = d3.scale.linear()
            .domain([0, data.length])
            .range([0, chart_width]);

        var y = d3.scale.linear()
            .domain([0, d3.max(data, function(d) { return d[metric]; })])
            .range([quality_chart_height, 0]);

        var xAxis = d3.svg.axis()
            .scale(x)
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient("left");

        var svg = d3.select('#region_coverage');

        if (svg.selectAll('rect').length == 0 || svg.selectAll('rect')[0].length == 0) {
            svg = d3.select('#region_coverage').append("svg")
            .attr("width", chart_width  + region_chart_margin.left + region_chart_margin.right)
            .attr("height", quality_chart_height + region_chart_margin.top + region_chart_margin.bottom)
            .append("g")
            .attr('id', 'inner_graph')
            .attr("transform", "translate(" + region_chart_margin.left + "," + region_chart_margin.top + ")");

            var bar = svg.selectAll(".bar")
                .data(data)
                .enter().append("g")
                .attr("class", "bar");

            bar.append("rect")
                .attr("x", function(d, i) { return x(i); })
                .attr("width", chart_width/data.length - 1)
                .attr("height", function(d) { return quality_chart_height - y(d[metric]); })
                .attr("y", function(d) { return y(d[metric]); });

            xAxis = d3.svg.axis()
                .scale(x)
                .tickFormat(function(d) { return ref[d]; })
                .innerTickSize(0)
                .orient("bottom");

            svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + quality_chart_height + ")")
                .call(xAxis);

            svg.append("g")
                .attr("class", "y axis")
                .call(yAxis);
        } else {
            svg = d3.select('#region_coverage').select('svg').select('#inner_graph');
            svg.select(".y.axis")
                .transition()
                .call(yAxis);

            svg.selectAll('rect')
                .data(data)
                .transition()
                .duration(500)
                .attr("x", function(d, i) { return x(i); })
                .attr("width", chart_width/data.length - 1)
                .attr("height", function(d) { return quality_chart_height - y(d[metric]); })
                .attr("y", function(d) { return y(d[metric]); });
        }
    } else {
        var data = {};
        $.each(raw_data[0], function(d, i) {
            var num = parseInt(d);
            if (!isNaN(num)) {
                data[d] = raw_data[0][d];
            }
        });

        var coverages = Object.keys(data);
        var all_labels = coverages;

        var chart_width = region_chart_width;
        var x = d3.scale.linear()
            .domain([0, coverages.length])
            .range([0, chart_width]);

        var y = d3.scale.linear()
            .domain([0, d3.max(coverages, function(d) { return data[d]; })])
            .range([quality_chart_height, 0]);

        var xAxis = d3.svg.axis()
            .scale(x)
            .tickFormat(function(d) { return all_labels[d - 1]; })
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient("left");

        svg = d3.select('#region_coverage').append("svg")
            .attr('id', 'inner_svg')
            .attr("width", chart_width + region_chart_margin.left + region_chart_margin.right)
            .attr("height", quality_chart_height + region_chart_margin.top + region_chart_margin.bottom)
            .append("g")
            .attr('id', 'inner_graph')
            .attr("transform", "translate(" + region_chart_margin.left + "," + region_chart_margin.top + ")");

        var bar = svg.selectAll(".bar")
            .data(coverages)
            .enter().append("g")
            .attr("class", "bar");

        bar.append("rect")
            .attr("x", function(d, i) { return x(i); })
            .attr("width", chart_width/coverages.length)
            .attr("height", function(d) { return quality_chart_height - y(data[d]); })
            .attr("y", function(d) { return y(data[d]); });

        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + quality_chart_height + ")")
            .call(xAxis)
            .selectAll("text")
            .attr("transform", "translate(0, 10) rotate(45)");

        svg.append("g")
            .attr("class", "y axis")
            .call(yAxis);

        svg.append("text")
            .attr("class", "x label")
            .attr("text-anchor", "middle")
            .style("font-size", "12px")
            .attr("x", region_chart_width/3)
            .attr("y", quality_chart_height + 50)
            .text(">= Coverage");
        svg.append("text")
            .attr("class", "y label")
            .attr("text-anchor", "middle")
            .style("font-size", "12px")
            .attr("transform", "rotate(-90)")
            .attr("x", -quality_chart_height/2)
            .attr("y", -40)
            .text("Fraction individuals covered");
    }
}


function add_line_to_quality_histogram(data, position, container, log) {
    //Takes dataset (for range) and datapoint and draws line in container
    //If line is already in container, transitions to new line
    var low_value = d3.min(data, function (d) { return d[0]; });
    var high_value = d3.max(data, function (d) { return d[0]; });
    if (log) {
        xscale = d3.scale.log()
            .domain([low_value, high_value])
            .range([0, quality_chart_width]);
    } else {
        xscale = d3.scale.linear()
            .domain([low_value, high_value])
            .range([0, quality_chart_width]);
    }
    x = function(d) {
        var pos;
        if (d > high_value) {
            pos = xscale(high_value);
        } else if (d < low_value) {
            pos = xscale(low_value);
        } else {
            pos = xscale(d);
        }
        return pos;
    };
    var svg = d3.select(container).select('svg').select('#inner_graph');
    if (svg.selectAll('.line').length == 0 || svg.selectAll('.line')[0].length == 0) {
        var lines = svg.selectAll(".line")
                    .data([position])
                    .enter().append("g")
                    .attr("class", "line");
        lines.append('line')
                .attr("x1", function(d) { return x(d); })
                .attr("x2", function(d) { return x(d); })
                .attr("y1", quality_chart_height)
                .attr("y2", 0)
                .attr("stroke-width", 2)
                .attr("stroke", "red");
    } else {
        svg.selectAll('.line').select('line')
            .data([position])
            .transition()
            .duration(500)
            .attr("x1", function(d) { return x(d); })
            .attr("x2", function(d) { return x(d); })
            .attr("y1", quality_chart_height)
            .attr("y2", 0)
            .attr("stroke-width", 2)
            .attr("stroke", "red");
    }
}

function draw_quality_histogram(data, container, log, integer_scale, xlabel, ylabel) {
    //Takes histogram data as a list of [midpoint, value] and puts into container
    //If data already in container, transitions to new data
    var x;
    if (log) {
        x = d3.scale.log();
    } else {
        x = d3.scale.linear();
    }
    x.domain([d3.min(data, function (d) {
        return d[0];
    }), d3.max(data, function (d) {
        return d[0];
    })])
        .range([0, quality_chart_width]);
    var bar_width = x(data[1][0]) - x(data[0][0]);
    var y = d3.scale.linear()
        .domain([d3.min(data, function(d) { return d[1]; }), d3.max(data, function(d) { return d[1]; })])
        .range([quality_chart_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom")
        .ticks(7);

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");
    if (integer_scale) {
        yAxis.tickFormat(d3.format("d"));
    }

    var svg = d3.select(container);
    if (svg.selectAll('rect').length == 0 || svg.selectAll('rect')[0].length == 0) {
        svg = d3.select(container).append("svg")
            .attr("width", quality_chart_width + quality_chart_margin.left + quality_chart_margin.right)
            .attr("height", quality_chart_height + quality_chart_margin.top + quality_chart_margin.bottom)
            .append("g")
            .attr('id', 'inner_graph')
            .attr("transform", "translate(" + quality_chart_margin.left + "," + quality_chart_margin.top + ")");
        svg.append("text")
            .attr("class", "x label")
            .attr("text-anchor", "middle")
            .style("font-size", "12px")
            .attr("x", quality_chart_width/2)
            .attr("y", quality_chart_height + xoffset)
            .text(xlabel);
        svg.append("text")
            .attr("class", "y label")
            .attr("text-anchor", "middle")
            .style("font-size", "12px")
            .attr("transform", "rotate(-90)")
            .attr("x", -quality_chart_height/2)
            .attr("y", -yoffset)
            .text(ylabel);
        var bar = svg.selectAll(".bar")
            .data(data)
            .enter().append("g")
            .attr("class", "bar");

        bar.append("rect")
            .attr("x", function(d) { return x(d[0]); })
            .attr("width", bar_width)
            .attr("height", function(d) { return quality_chart_height - y(d[1]); })
            .attr("y", function(d) { return y(d[1]); });

        if (container == '#quality_metric_container') {
            svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + quality_chart_height + ")")
                .style("font-size", "10px")
                .call(xAxis)
                .selectAll("text")
                .style("text-anchor", "end")
                .attr("dx", "-.8em")
                .attr("dy", ".15em")
                .attr("transform", function(d) {
                    return "rotate(-45)"
                });
            svg.append("g")
                .attr("class", "y axis")
                .style("font-size", "10px")
                .call(yAxis);
        } else {
            svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + quality_chart_height + ")")
                .call(xAxis);
            svg.append("g")
                .attr("class", "y axis")
                .call(yAxis);
        }
    } else {
        svg = d3.select(container).select('svg').select('#inner_graph');

        if (container == '#quality_metric_container') {
            svg.select(".x.axis")
                .transition()
                .attr("transform", "translate(0," + quality_chart_height + ")")
                .call(xAxis)
                .selectAll("text")
                .style("text-anchor", "end")
                .attr("dx", "-.8em")
                .attr("dy", ".15em")
                .attr("transform", function(d) {
                    return "rotate(-45)"
                });
        } else {
            svg.select(".x.axis")
            .transition()
            .attr("transform", "translate(0," + quality_chart_height + ")")
            .call(xAxis);
        }

        svg.select(".y.axis")
            .transition()
            .call(yAxis);

        svg.select('.x.label')
            .text(xlabel);
        svg.select('.y.label')
            .text(ylabel);
        svg.selectAll('rect')
            .data(data)
            .transition()
            .duration(500)
            .attr("x", function(d) { return x(d[0]); })
            .attr("width", bar_width)
            .attr("height", function(d) { return quality_chart_height - y(d[1]); })
            .attr("y", function(d) { return y(d[1]); });
    }
}


$(document).ready(function() {
    draw_region_coverage(window.base_coverage, 'mean', window.variant.ref);
    $('.coverage_metric_buttons').change(function () {
        var v = $(this).attr('id').replace('_covmet_button', '');
        $('.coverage_subcat_selectors').hide();
        if (v == 'covered') {
            $('#over_x_select_container').show();
            v = $('#over_x_select').val();
        } else {
            $('#average_select_container').show();
            v = $("#average_select").val();
        }
        draw_region_coverage(window.base_coverage, v, window.variant.ref);
    });
    $('#over_x_select').change(function () {
        draw_region_coverage(window.base_coverage, $(this).val(), window.variant.ref);
    });
    $('#average_select').change(function () {
        draw_region_coverage(window.base_coverage, $(this).val(), window.variant.ref);
    });
});

function check_for_variant_in_clinvar() {
    var clinvar_searches = _.map(variant.rsids, function(rsid) {
        var clinvar_query = 'term={RSID}[Variant%20ID]&retmode=json'.replace('{RSID}', rsid);
        return {
            xhr_url: 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&' + clinvar_query,
            name: rsid,
            webpage_url: 'http://www.ncbi.nlm.nih.gov/clinvar?' + clinvar_query
        };
    });
    var clinvar_query = 'term={CHROM}[Chromosome]%20AND%20{POS}[Base%20Position%20for%20Assembly%20GRCh37]&retmode=json'
        .replace('{CHROM}', variant.chrom).replace('{POS}', variant.pos);
    clinvar_searches.push({
        xhr_url: 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&' + clinvar_query,
        name: variant.chrom + ':' + variant.pos,
        webpage_url: 'http://www.ncbi.nlm.nih.gov/clinvar?' + clinvar_query
    });

    _.each(clinvar_searches, function(clinvar_search) {
        $.getJSON(clinvar_search.xhr_url)
        .done(function(data) {
            var link_text = (data.esearchresult.count !== "0") ? 'Open <%= name %> in ClinVar' : '<%= name %> is not in ClinVar';
            var link_style = 'style="float:left; clear:both"';
            if ($('#clinvar_loading').length !== 0) {
                $('#clinvar_loading').remove();
                link_style = ''; // One of the links must not float or clear, so that it can sit next to <dt>ClinVar</dt>
            }
            $('#clinvar').append(_.template(
                '<a href="<%= url %>" target="_blank" <%= style %>>' + link_text + ' <i class="fa fa-external-link"></i> </a> '
            )({name: clinvar_search.name, url: clinvar_search.webpage_url, style: link_style}));
        });
    });
}
check_for_variant_in_clinvar();

$(function() {
    $('.transcript_toggle').on('click', function(e) {
        e.preventDefault();
        var $collapse = $(this).closest('.collapse-group').find('.collapse');
        $collapse.collapse('toggle');

        //change text
        $(this).find('.transcript-toggle-text.hidden').removeClass('hidden').hide();
        $(this).find('.transcript-toggle-text').each(function() { $(this).toggle(); });
    });
});

var af_buckets = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
function get_af_bucket_text(bin) {
    if (bin == 'singleton' || bin == 'doubleton') {
        return 'This is the site quality distribution for all ' + bin + 's in ' + window.dataset_name + '.';
    } else if (bin == '0.0001') {
        return 'This is the site quality distribution for all variants with AF < ' + bin + ' in ' + window.dataset_name + '.';
    } else {
        return 'This is the site quality distribution for all variants with ' + af_buckets[af_buckets.indexOf(parseFloat(bin)) - 1] + ' < AF < ' + bin + ' in ' + window.dataset_name + '.';
    }
}

$(document).ready(function() {
    $('#pop_afs_table').tablesorter({
        sortList: [[0,0]],
    });

    if (window.variant != null && 'genotype_depths' in window.variant) {
        draw_quality_histogram(window.variant.genotype_depths[0], '#quality_display_container', false, true, 'Depth', 'Number of Individuals');

        $('.quality_display_buttons, .quality_full_site_buttons').change(function() {
            setTimeout(function() {
                var v = $('.quality_display_buttons.active').attr('id').replace('_button', '');
                var f = $('.quality_full_site_buttons.active').attr('id') == 'variant_site_button' ? 0 : 1;
                var ylab = 'Number of ' + (f ? 'Variant Carriers' : 'Individuals');
                var xlab = $('.quality_display_buttons.active').text().replace(/^\s+|\s+$/g, '');
                draw_quality_histogram(window.variant[v][f], '#quality_display_container', false, xlab === 'Depth', xlab, ylab);
            }, 0);
        });

        // Quality metric histograms
        var data = _.zip(_.map(window.metrics['Site Quality']['mids'], Math.exp), window.metrics['Site Quality']['hist']);
        draw_quality_histogram(data, '#quality_metric_container', true, false, 'Site Quality', 'Variants');
        var bin = window.metrics['Site Quality']['metric'].split('_')[1];
        $('#site_quality_note').html(get_af_bucket_text(bin));
        var pos = $('#quality_metric_select').val().split(': ')[1];
        add_line_to_quality_histogram(data, pos, '#quality_metric_container', true);
        var log_scale_metrics = ['Site Quality', 'Total Depth'];
        $('#quality_metric_select').change(function() {
            var v = $(this).val().split(': ');
            var log = false;
            $('#site_quality_note').html('');
            var data;
            if (log_scale_metrics.indexOf(v[0]) > -1) {
                data = _.zip(_.map(window.metrics[v[0]]['mids'], Math.exp), window.metrics[v[0]]['hist']);
                log = true;
            } else {
                data = _.zip(window.metrics[v[0]]['mids'], window.metrics[v[0]]['hist']);
            }
            if (v[0] == 'Site Quality') {
                var bin = window.metrics['Site Quality']['metric'].split('_')[1];
                $('#site_quality_note').html(get_af_bucket_text(bin));
            }
            draw_quality_histogram(data, '#quality_metric_container', log, false, v[0], 'Variants');
            add_line_to_quality_histogram(data, v[1], '#quality_metric_container', log);
        });
    } else {
        $('#quality_metrics_container').hide();
    }
});
