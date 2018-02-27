function draw_region_coverage(raw_data, metric, ref) {
    // pjvh removed a bunch of functionality from this.  There's some useless code left behind.
    // If this function gets a single base, it draws the full distribution.
    // If it receives multiple bases, it draws a coverage graph, letting the user select mean, median, % > X
    // TODO: always draw coverage for only the first base.
    region_chart_width = 500;
    region_chart_margin = {top: 10, right: 50, bottom: 55, left: 70};
    if (raw_data.length > 1) {
        var data = raw_data;
        var chart_width = _.min([region_chart_width, data.length*30]);
        //var x = d3.scale.linear()
        var x = d3.scaleLinear()
            .domain([0, data.length])
            .range([0, chart_width]);

        //var y = d3.scale.linear()
        var y = d3.scaleLinear()
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
                .call(yAxis);

            svg.selectAll('rect')
                .data(data)
                .attr("x", function(d, i) { return x(i); })
                .attr("width", chart_width/data.length - 1)
                .attr("height", function(d) { return quality_chart_height - y(d[metric]); })
                .attr("y", function(d) { return y(d[metric]); });
        }
    } else if (raw_data.length == 1) {
        var data = {};
        $.each(raw_data[0], function(d, i) {
            var num = parse_int(d);
            if (!isNaN(num)) {
                data[d] = raw_data[0][d];
            }
        });

        var coverages = Object.keys(data);
        var all_labels = coverages;

        var chart_width = region_chart_width;
        //var x = d3.scale.linear()
        var x = d3.scaleLinear()
            .domain([0, coverages.length])
            .range([0, chart_width]);

        //var y = d3.scale.linear()
        var y = d3.scaleLinear()
            .domain([0, d3.max(coverages, function(d) { return data[d]; })])
            .range([quality_chart_height, 0]);


        //var xAxis = d3.svg.axis()
        var xAxis = d3.axisBottom()
            .scale(x)
            .tickFormat(function(d) { return all_labels[d - 1]; });
            //.orient("bottom");

        //var yAxis = d3.svg.axis()
        var yAxis = d3.axisLeft()
            .scale(y)
            //.orient("left")
            .ticks(5)
            .tickFormat(d3.format('%'));

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


function draw_metric_histograms(dom_element, histogram1_bins, histogram2_bins, bin_width, label1, label2, x_label, y_label, title) {
    // parameters
    var svg_width = 250;
    var svg_height = 280;
    var font_size = 12;
    var margin = { top: 40, right: 10, bottom: 35, left: 55 };

    const color1 = "#4169E1";
    const color2 = "#ff69b4";

    /* create SVG element */
    var svg = d3.select(dom_element.get(0)).append("svg")
        .attr("width", svg_width)
        .attr("height",svg_height);

    /* create X label */
    var x_label_text = svg.append("text")
        .attr("x", svg_width / 2)
        .attr("y", 0)
        .style("text-anchor", "middle")
        .style("dominant-baseline", "hanging")
        .attr("font-size", font_size)
        .text(x_label);    
    x_label_text.attr("y", svg_height - x_label_text.node().getBBox().height);

    /* create Y label */
    var y_label_text = svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("x", - (svg_height / 2))
        .attr("y", 0)
        .style("text-anchor", "middle")
        .style("dominant-baseline", "hanging")
        .attr("font-size", font_size)
        .text(y_label);

    /* create title */
    if (title) {
        var title_text = svg.append("text")
            .attr("x", svg_width / 2)
            .attr("y", 0)
            .style("text-anchor", "middle")
            .style("dominant-baseline", "hanging")
            .attr("font-size", font_size * 1.5)
            .text(title);
    }

    /* Draw legend */
    function add_legend_label(svg_element, x, y, label, color) {
        svg_element.append("rect")
            .attr("x", x)
            .attr("y", y)
            .attr("width", font_size + 1)
            .attr("height", font_size + 1)
            .attr("fill", color)
            .attr("fill-opacity", "0.50");
        return y + svg_element.append("text")
            .attr("x", font_size + 3)
            .attr("y", y)
            .attr("text-anchor", "start")
            .attr("dominant-baseline", "text-before-edge")
            .attr("font-size", font_size)
            .text(label).node().getBBox().height;
    }

    var legend = svg.append("g");
    if (title) {
        legend.attr("transform", "translate(" + 0 + "," + title_text.node().getBBox().height + ")");
    }
    add_legend_label(legend, 0, add_legend_label(legend, 0, 0, label1, color1), label2, color2);

    /* Set width and height for the plotting area */
    var width = svg_width - margin.left - margin.right;
    var height = svg_height - margin.top - margin.bottom;

    /* Transform counts to proportion */

    const total1 = d3.sum(histogram1_bins, function(d) { return d[1]; });
    const total2 = d3.sum(histogram2_bins, function(d) { return d[1]; });
    const y_extent1 = d3.extent(histogram1_bins, function(d) { return d[1] / total1; });
    const x_extent1 = d3.extent(histogram1_bins, function(d) { return d[0]; });
    const y_extent2 = d3.extent(histogram2_bins, function(d) { return d[1] / total2; });
    const x_extent2 = d3.extent(histogram2_bins, function(d) { return d[0]; });
    
    /* Set X and Y scales */
    var scale_x = d3.scaleLinear();
    scale_x.domain([(x_extent2[0] > x_extent1[0] ? x_extent2[0] : x_extent1[0]) - bin_width / 2, (x_extent2[1] > x_extent1[1] ? x_extent2[1] : x_extent1[1])  + bin_width / 2]);
    scale_x.range([0, width]);

    var scale_y = d3.scaleLinear();
    scale_y.domain([0, y_extent2[1] > y_extent1[1] ? y_extent2[1] : y_extent1[1]]);
    scale_y.range([height, 0]);

    svg = svg.append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    svg.append("g").selectAll("rect")
        .data(histogram1_bins)
        .enter()
        .append("rect")
        .attr("x", function(d) { return scale_x(d[0]) - scale_x(bin_width / 2); })
        .attr("y", function(d) { return scale_y(d[1] / total1); })
        .attr("width", scale_x(bin_width))
        .attr("height", function(d) { return scale_y(0) - scale_y(d[1] / total1); })
        .attr("shape-rendering", "crispEdges") // to avoid white space between adjacent rectangles
        .attr("stroke", "none")
        .attr("stroke-width", "0")
        .attr("fill", color1)
        .attr("fill-opacity", "0.50")

    svg.append("g").selectAll("rect")
        .data(histogram2_bins)
        .enter()
        .append("rect")
        .attr("x", function(d) { return scale_x(d[0]) - scale_x(bin_width / 2); })
        .attr("y", function(d) { return scale_y(d[1] / total2); })
        .attr("width", scale_x(bin_width))
        .attr("height", function(d) { return scale_y(0) - scale_y(d[1] / total2); })
        .attr("shape-rendering", "crispEdges") // to avoid white space between adjacent rectangles
        .attr("stroke", "none")
        .attr("stroke-width", "0")
        .attr("fill", color2)
        .attr("fill-opacity", "0.50");

    /* Set up axes */
    var x_axis = d3.axisBottom(scale_x);
    svg.append("g")
        .attr("class", "x-axis")
        .attr("transform", "translate(0," + height + ")")
        .call(x_axis);

    var y_axis = d3.axisLeft(scale_y);
    svg.append("g")
        .attr("class", "y-axis")
        .call(y_axis);


    /* Setup tooltip */
    var tooltip = d3.tip()
        .attr("class", "d3-tip")
        .offset([-8, 0])
        .html(function(d) {
            var html = "";
            html += "All Individuals: " + d[0][1];
            html += "<br>Variant Carriars: " + d[1][1];
            return html; 
        });
    svg.call(tooltip);

    /* Create focus line */
    const focus = svg.append("g")
        .style("display", "none");

    focus.append("path")
        .attr("d", function(d) { return " M 0 0 L 0 " + height;  })
        .style("stroke", "red")
        .style("stroke-width", "1");


    /* Draw transparent rectangle to catch mouse events */
    svg.append("rect")
        .attr("width", width)
        .attr("height", height)
        .style("fill", "none")
        .style("pointer-events", "all")
        .on('mouseover', function() { focus.style("display", null);})
        .on('mouseout', function() { focus.style("display", "none"); tooltip.hide();})
        .on('mousemove', mousemove );

    const bisect = d3.bisector(function(d) { return d[0]; }).left;

    function closest_bin(bins, mouse_x) {
        const i = bisect(bins, mouse_x, 1);
        const xy0 = bins[i - 1];
        if (i < bins.length) {
            const xy1 = bins[i];
            return mouse_x - xy0[0] > xy1[0] - mouse_x ? xy1 : xy0;
        } else {
            return xy0;
        }
    }

    function mousemove() {
        const x0 = scale_x.invert(d3.mouse(this)[0]);
        bin1 = closest_bin(histogram1_bins, x0);
        bin2 = closest_bin(histogram2_bins, x0);
        focus.attr('transform', "translate(" + scale_x(x0) + ", 0)");
        tooltip.show([bin1, bin2], focus.node());
    }

}

$(document).ready(function() {
    /* draw_region_coverage(window.base_coverage, 'mean', window.variant.ref);
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
    });*/
});

function check_for_variant_in_clinvar() {
    var clinvar_searches = _.map(variant.rsids, function(rsid) {
        var clinvar_query = 'term={RSID}&retmode=json'.replace('{RSID}', rsid);
        return {
            xhr_url: 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&' + clinvar_query,
            name: rsid,
            webpage_url: 'https://www.ncbi.nlm.nih.gov/clinvar?' + clinvar_query
        };
    });
    var clinvar_query = 'term={CHROM}[Chromosome]%20AND%20{POS}[Base%20Position%20for%20Assembly%20GRCh37]&retmode=json'
        .replace('{CHROM}', variant.chrom).replace('{POS}', variant.pos);
    clinvar_searches.push({
        xhr_url: 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&' + clinvar_query,
        name: variant.chrom + ':' + variant.pos,
        webpage_url: 'https://www.ncbi.nlm.nih.gov/clinvar?' + clinvar_query
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
        return 'This is the site quality distribution for all ' + bin + 's in ' + window.model.dataset_name + '.';
    } else if (bin == '0.0001') {
        return 'This is the site quality distribution for all variants with AF < ' + bin + ' in ' + window.model.dataset_name + '.';
    } else {
        return 'This is the site quality distribution for all variants with ' + af_buckets[af_buckets.indexOf(parseFloat(bin)) - 1] + ' < AF < ' + bin + ' in ' + window.model.dataset_name + '.';
    }
}

function create_pop_afs_table() {
    $('#pop_afs_table').DataTable({
        paging: false,
        searching: false,
        info: false,
        ordering: false,
    });
}

function draw_percentiles_legend(dom_element) {
    // parameters
    var svg_width = 180;
    var svg_height = 45;
    var font_size = 12;
    var margin = { top: 15, right: 8, bottom: 18, left: 7 };

    /* create SVG element */
    var svg = d3.select(dom_element.get(0)).append("svg")
        .attr("width", svg_width)
        .attr("height", svg_height);

    /* create title */
    var text = svg.append("text")
        .attr("x", svg_width / 2)
        .attr("y", 0)
        .style("text-anchor", "middle")
        .style("dominant-baseline", "text-before-edge")
        .attr("font-size", font_size)
        .text("% of PASS variants");

    /* Set width and height for the plotting area */
    var width = svg_width - margin.left - margin.right;
    var height = svg_height - margin.top - margin.bottom;

    svg = svg.append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    /* Set X scale */
    var scale_x = null; 
    scale_x = d3.scaleLinear();
    scale_x.domain([0, 100]);
    scale_x.range([0, width - 1]);

    /* Set color scale */
    var colorScale = d3.scaleSequential(d3.interpolateRdYlBu).domain([0, 100]);

    svg.selectAll("rect")
        .data(d3.range(width))
        .enter()
        .append("rect")
        .attr("x", function(d, i) { return i; })
        .attr("y", 0)
        .attr("height", height)
        .attr("width", 1)
        .attr("shape-rendering", "crispEdges")
        .style("fill", function(d, i) { return colorScale(d); });

    /* Draw bounding box */
    svg.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", width)
        .attr("height", height)
        .style("pointer-events", "none") // to prevent bubbling
        .style("stroke", "#434343")
        .style("stroke-width", "1")
        .style("fill", "none");

    /* Set up axes */
    var x_axis = d3.axisBottom(scale_x).ticks(5);
    svg.append("g")
        .attr("class", "x-axis")
        .attr("transform", "translate(0," + height + ")")
        .call(x_axis);
}

function draw_metric_percentiles(dom_element, percentiles, value, value_p_low, value_p_high) {
    // parameters
    var svg_width = 180;
    var svg_height = 30;

    //var font_size = 12;
    var margin = { top: 2, right: 7, bottom: 8, left: 7 };

    /* create SVG element */
    var svg = d3.select(dom_element.get(0)).append("svg")
        .attr("width", svg_width)
        .attr("height", svg_height);

    /* Set width and height for the plotting area */
    var width = svg_width - margin.left - margin.right;
    var height = svg_height - margin.top - margin.bottom;

    svg = svg.append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    /* Set X scale */
    var scale_x = null; 
    scale_x = d3.scaleLinear();
    scale_x.domain([0, d3.max(percentiles, function(d) { return d.probability; })]);
    scale_x.range([0, width]);

    /* Prepare data for easier visualization */
    percentiles.sort(function(o1, o2) { return o1.probability < o2.probability ? -1 : 1 }); // sort percentile by probability
    var previous_probability = 0;
    var previous_value = null;
    var previous_n = 0;
    var previous_n_pass = 0;
    percentiles.forEach(function(percentile) {
        percentile.previous_probability = previous_probability;
        percentile.previous_value = previous_value;
        percentile.bin_n = percentile.n - previous_n;
        percentile.bin_n_pass = percentile.n_pass - previous_n_pass;
        previous_probability = percentile.probability;
        previous_value = percentile.value;
        previous_n = percentile.n;
        previous_n_pass = percentile.n_pass;
    });

    /* Setup tooltip */
    var tooltip = d3.tip()
        .attr("class", "d3-tip")
        .offset([-8, 0])
        .html(function(d) {
            var html = ""
            html += "From <i>P</i><sub>" + (d.previous_probability * 100)  + "</sub> = " + d.previous_value;
            html += "<br>To <i>P</i><sub>" + (d.probability * 100) + "</sub> = " + d.value;
            html += "<br>PASS (%) = " + (d.bin_n_pass / d.bin_n * 100).toFixed(4);
            return html; 
        });
    svg.call(tooltip);

    /* Draw percentiles */
    svg.selectAll("rect")
        .data(percentiles)
        .enter()
        .append("rect")
        .attr("x", function(d) { return scale_x(d.previous_probability); })
        .attr("y", 0)
        .attr("width", function(d) { return scale_x(d.probability) - scale_x(d.previous_probability); })
        .attr("height", height)
        .style("stroke", "none")
        .style("stroke-width", "0")
        .attr("shape-rendering", "crispEdges") // to avoid white space between adjacent rectangles
        .style("fill", function(d) {
            if (d.bin_n > 0) {
                return d3.interpolateRdYlBu(d.bin_n_pass / d.bin_n);
            } else {
                return "#a1a1a1";
            }
        })
        .on("mouseover", function(d) {
            focus.style("display", null)
            focus.select("rect")
                .attr("x", scale_x(d.previous_probability))
                .attr("width", scale_x(d.probability) - scale_x(d.previous_probability));
            tooltip.show(d);
        })
        .on("mouseout", function(d) {
            focus.style("display", "none");
            tooltip.hide();
        });

    /* Draw percentile rectangle mark if difference between percentile values is more than 2px. Otherwise, draw triangle. */
    if (Math.round(scale_x(value_p_high)) - Math.round(scale_x(value_p_low)) >= 3) { 
        const p_rect = svg.append("g")
            .attr("transform", function(d) { return "translate(" + scale_x(value_p_low)  + "," + height + ")"; });
        p_rect.append("path")
            .attr("d", function(d) {
                const length = (scale_x(value_p_high) - scale_x(value_p_low));
                return " M 0 0 L 0 7 L " + length + " 7 L " + length + " 0 L 0 0";
            })
            .style("stroke", "none")
            .style("stroke-width", "0")
            .attr("shape-rendering", "crispEdges") // to avoid white space between adjacent triangle
            .style("fill", "red");
    } else {
        const p_tri = svg.append("g")
            .attr("transform", function(d) { return "translate(" + scale_x(value_p_low)  + "," + height + ")"; });
  
        p_tri.append("path")
           .attr("d", function(d) {
               return " M 0 0 L -7 7 L 7 7 L 0 0";
           })
           .style("stroke", "none")
           .style("stroke-width", "0")
           .attr("shape-rendering", "crispEdges") // to avoid white space between adjacent triangle
           .style("fill", "red");
    }

    /* Setup focus rectangle to highlight selected percentile */
    const focus = svg.append("g")
        .style("display", "none");

    focus.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", 0)
        .attr("height", height)
        .attr("shape-rendering", "crispEdges") // to avoid white space between adjacent rectangles
        .style("pointer-events", "none") // to prevent bubbling
        .attr("fill", "#000000")
        .attr("fill-opacity", "0.5");

    /* Draw bounding box */
    svg.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", width)
        .attr("height", height)
        .style("pointer-events", "none") // to prevent bubbling
        .style("stroke", "#434343")
        .style("stroke-width", "1")
        .style("fill", "none");
}


$(document).ready(function() {
    create_pop_afs_table();

      $('#variant-nav').affix({ 
          offset: {
              top: function() { return $('#variant-name').height() + $('#variant-name').offset().top; }
          }
      });


    if (window.variant != null && 'genotype_depths' in window.variant) {
        var dom_element = null;
        
        dom_element = $("<div></div>").appendTo("#sequence-depth-plots");
        dom_element.addClass("col-lg-3 col-md-4 col-sm-6 col-xs-12");
        draw_metric_histograms(dom_element, window.variant.genotype_depths[0], window.variant.genotype_depths[1], 5, "All Individuals", "Variant Carriers", "Sequence Depth", "Proportion of Individuals", "");

        dom_element = $("<div></div>").appendTo("#genotype-quality-plots");
        dom_element.addClass("col-lg-3 col-md-4 col-sm-6 col-xs-12");
        draw_metric_histograms(dom_element, window.variant.genotype_qualities[0], window.variant.genotype_qualities[1], 5, "All Individuals", "Variant Carriers", "Genotype Quality", "Number of Individuals", "");

        draw_percentiles_legend($("#site-metrics-legend"));
   
        window.metrics.forEach(function (metric) {
            var percentiles = null;
            var value = null;
            if ((metric.name in window.variant.quality_metrics_percentiles) && (metric.name in window.variant.quality_metrics)) {
                percentiles =  window.variant.quality_metrics_percentiles[metric.name];
                value = window.variant.quality_metrics[metric.name];
                dom_element = $("<tr></tr>").appendTo("#site-metrics-plots");
                $("<td></td>").appendTo(dom_element).html(metric.description);
                $("<td></td>").appendTo(dom_element).html(value);
                $("<td></td>").appendTo(dom_element).html((percentiles[1] * 100).toFixed(2));
                draw_metric_percentiles($("<td></td>").appendTo(dom_element), metric.percentiles, value, percentiles[0], percentiles[1]);
            }
        });
    } else {
        $('#quality_metrics_container').hide();
    }
});
