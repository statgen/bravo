var total_width = $(window).width() < 768 ? $(window).width() : $(window).width() * 10 / 12;

var gene_chart_margin = {top: 10, right: 10, bottom: 5, left: 30};
if ($(window).width() < 768) {
    gene_chart_margin.left = 10;
}

var gene_chart_margin_lower = {
    top: 5,
    right: gene_chart_margin.right,
    bottom: 5,
    left: gene_chart_margin.left
};

var gene_chart_width = total_width - gene_chart_margin.left - gene_chart_margin.right;

var lower_gene_chart_height = 50 - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom;
var gene_chart_height = 300 - gene_chart_margin.top - gene_chart_margin.bottom - lower_gene_chart_height - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom;


// Padding can have variants and coverage.
// Margin is just blank.
var EXON_PADDING = 15;
var EXON_MARGIN = 60;

/*
    The following methods are for working with "Coding Coordinates",
    a coordinate space that we use to plot data in the coding regions of a transcript.

    Conceptually: suppose you lined up all the coding regions in a transcript, with some padding on each side,
    then plotted any variants that overlap. The position of a variant on this line is the coding position -
    this will obviously differ from the actual genomic or cds coordinates.

    Random notes:
        - coding coordinates have no concept of a gene - they are solely a property of a transcript
        - should probably have a map between coding coodinates and protein position
        - only "CDS" is shown by default.  "UTR" and "exon" are hidden until "show UTRs" is checked.
 */
window.get_position_mapping = _.memoize(function(skip_utrs) {
    // Uses window.exons_and_utrs
    // Returns like [
    //    {real_start: 36649935, scaled_start: 0, length: 178},
    //    {real_start: 36650909, scaled_start: 179, length: 212}
    // ]

    var exons = window.exons_and_utrs;
    if (skip_utrs) {
        exons = _.filter(exons, function(exon) {
            return exon.feature_type === "CDS";
        });
        if (exons.length === 0) {
            exons = window.exons_and_utrs;
        }
    }
    if (exons.length === 0) {
        return [];
    }

    // In theory, overlap is fine, so long as both map from the same real coordinates to scaled coordinates.
    // That is, as long as the difference between read_start is the same as the difference between scaled_start.
    // In practice, overlap brings out the bugs that I haven't fixed, so I'll unify overlapping sections.
    exons = _.sortBy(exons, "start");
    var pos_mapping = [{
        real_start: exons[0].start - EXON_PADDING,
        scaled_start: EXON_MARGIN,
        length: exons[0].stop - exons[0].start + EXON_PADDING*2,
    }];
    for (var i=1; i<exons.length; i++) {
        var prev_map = pos_mapping[pos_mapping.length-1];
        var new_map = {
            real_start: exons[i].start-EXON_PADDING,
            length: exons[i].stop - exons[i].start + EXON_PADDING*2
        };
        var real_end_of_prev_map = prev_map.real_start + prev_map.length;
        var scaled_end_of_prev_map = prev_map.scaled_start + prev_map.length + 1; // +1 ?
        var scaled_gap_between_maps = Math.min(new_map.real_start - real_end_of_prev_map, EXON_MARGIN*2);
        new_map.scaled_start = scaled_end_of_prev_map + scaled_gap_between_maps;

        if (scaled_gap_between_maps <= 0) {
            var possible_new_length = new_map.scaled_start + new_map.length - prev_map.scaled_start;
            prev_map.length = Math.max(prev_map.length, possible_new_length);
        } else {
            pos_mapping.push(new_map);
        }
    }
    return pos_mapping;
});

window.get_coding_coordinate = function(position, skip_utrs) {
    var pos_mapping = window.get_position_mapping(skip_utrs);
    // TODO: binary search
    for (var i=0; i<pos_mapping.length; i++) {
        var m = pos_mapping[i];
        if (position < m.real_start) {
            return null;
        } else if (position <= m.real_start + m.length) {
            return m.scaled_start + position - m.real_start;
        }
    }
    return null;
};

window.precalc_coding_coordinates_for_bin = function(bin, skip_utrs) {
    // TODO: if a bin extends across two exons, then it needs to be split.
    //       to split a bin, this function should be `window.get_bins_with_coding_coordinates = function(bins, skip_utrs)`
    //       and it should `rv=[]`, `rv.append`, etc.

    // bin is like {start:2162705, end:2162706, ...}
    // return is like {start:2162705, start_coding_noutr:50, start_coding:50, end:2162706, end_coding_noutr:51, end_coding:51, ...}

    var pos_mapping = window.get_position_mapping(skip_utrs);

    var key_suffix = skip_utrs ? '_coding_noutr' : '_coding';
    for (var i=0; i<pos_mapping.length; i++) {
        var m = pos_mapping[i];
        if (bin.end < m.real_start) {
            bin['start'+key_suffix] = null;
            bin['end'+key_suffix] = null;
            return;
        } else if (bin.start <= m.real_start + m.length) {
            // +1 and -1 just help make a little break
            bin['start'+key_suffix] = Math.max(m.scaled_start + 1, m.scaled_start + bin.start - m.real_start);
            bin['end'+key_suffix] = Math.min(m.scaled_start + m.length - 1, m.scaled_start + bin.end - m.real_start);
            return;
        }
    }
};

window.get_maximum_coding_position = _.memoize(function(skip_utrs) {
    var pos_mapping = window.get_position_mapping(skip_utrs);
    if (pos_mapping.length === 0) {
        return 0;
    }
    //assume that we start at 0 and go EXON_PADDING+EXON_MARGIN beyond the end of the last pos_mapping.
    //Note that nothing can actually reach that into the EXON_MARGIN, so this function is kinda lying.
    return pos_mapping[pos_mapping.length-1].length + pos_mapping[pos_mapping.length-1].scaled_start + EXON_PADDING + EXON_MARGIN;
});

window.precalc_coding_coordinates = function(objects) {
    _.each(objects, function(o) {
        o.pos_coding = get_coding_coordinate(o.pos, false)
        o.pos_coding_noutr = get_coding_coordinate(o.pos, true)
    });
};

window.precalc_coding_coordinates_for_bins = function(bins) {
    _.each(bins, function(bin) {
        window.precalc_coding_coordinates_for_bin(bin, false);
        window.precalc_coding_coordinates_for_bin(bin, true);
    });
};


function create_charts() {
    var skip_utrs = true;
    var coords = skip_utrs ? 'pos_coding_noutr' : 'pos_coding';
    var scale_type = 'overview';
    var maximum_coding_position = get_maximum_coding_position(skip_utrs);
    var chart_width = (scale_type == 'overview') ?
        gene_chart_width :
        Math.max(gene_chart_width, Math.min(gene_chart_width*4, maximum_coding_position*2));
    var exon_x_scale = d3.scale.linear()
        .domain([0, maximum_coding_position])
        .range([0, chart_width]);

    var coverage_track = d3.select('#gene_plot_container').append("svg")
        .attr("width", chart_width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .attr('id', 'coverage_svg')
        .attr('class', 'hidden-xs')
        .append("g")
        .attr('id', 'coverage_track')
        .attr("transform", "translate(" + gene_chart_margin.left + "," + gene_chart_margin.top + ")");
    d3.select('#gene_plot_container').append("br"); //make sure exon_track is below graph and not to the right of it
    var exon_track = d3.select('#gene_plot_container').append("svg")
        .attr("width", chart_width + gene_chart_margin_lower.left + gene_chart_margin_lower.right)
        .attr("height", lower_gene_chart_height)
        .attr('id', 'exon_svg')
        .append("g")
        .attr('id', 'exon_track')
        .attr("transform", "translate(" + gene_chart_margin_lower.left + "," + 0 + ")");

    // show variant category on hover
    window.variant_plot_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        if (d.category) {
            var csq = d.major_consequence.replace('_variant', '')
                    .replace('_', ' ')
                    .replace('utr', 'UTR')
                    .replace('3 prime', "3'")
                    .replace('5 prime', "5'")
                    .replace('nc ', "non-coding ");
            var output = csq + '<br/>' + d.chrom + ':' + d.pos.toLocaleString() + ' ' + d.ref + '&#8594;' + d.alt;
            if (d.major_consequence == 'missense_variant' || d.major_consequence == 'synonymous_variant') {
                output += '<br/>' + d.HGVSp;
            }
            output += '<br/>Frequency: ' + d.allele_freq.toPrecision(3);
            return output;
        } else {
            return 'None';
        }
    });
    coverage_track.call(window.variant_plot_tip);

    create_coverage_chart(coverage_track, coords, chart_width, exon_x_scale);
    create_exon_plot(exon_track, exon_x_scale, skip_utrs);
    create_variants_plot(exon_track, exon_x_scale, coords);
}

function create_coverage_chart(coverage_track, coords, chart_width, exon_x_scale) {
    var metric = 'mean';
    var coords_start = coords.replace('pos_', 'start_');
    var coords_end = coords.replace('pos_', 'end_');

    var area = d3.svg.area()
        .x( function(d) {
            return exon_x_scale(d['pos']);
        }).y0( function(d) {
            return gene_chart_height;
        }).y1( function(d) {
            return (metric in d) ? y(d[metric]) : gene_chart_height;
        });

    //This is a hack, revealing that the surrounding code needs changes.
    my_datum = [];
    _.each(window.coverage_stats, function(d) {
        if (d[coords_start] && d[coords_end]) {
            var o = [{}, {}, {}, {}];
            o[0]['pos'] = o[1]['pos'] = d[coords_start];
            o[2]['pos'] = o[3]['pos'] = 1 + d[coords_end];
            o[1][metric] = o[2][metric] = d[metric] || 0;
            o[0][metric] = o[3][metric] = 0;
            my_datum.push.apply(my_datum, o); //concatenate o to my_datum
        }
    });

    var max_cov = 1;
    if (metric == 'mean' || metric == 'median') {
        max_cov = d3.max(my_datum, function(d) { return d[metric]; });
    }
    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([gene_chart_height, 0]);

    coverage_track.append("path")
        .datum(my_datum)
        .style("fill", "steelblue")
        .attr('class', 'area')
        .attr("d", area);

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    coverage_track.append("g")
        .attr("class", "y axis")
        .call(yAxis);
}

function create_exon_plot(exon_track, exon_x_scale, skip_utrs) {
    var exon_color = "lightsteelblue";
    exon_track.selectAll("line.padded_exon")
        .data(window.exons_and_utrs)
        .enter()
        .append('line')
        .attr("class", "padded_exon")
        .attr("y1", lower_gene_chart_height/2)
        .attr("y2", lower_gene_chart_height/2)
        .attr("stroke-width", 5)
        .attr("stroke", exon_color);

    // plot exon rects
    exon_track.selectAll("rect.exon_or_utr")
        .data(window.exons_and_utrs)
        .enter()
        .append("rect")
        .attr('class', 'exon_or_utr')
        .style("fill", exon_color)
        .attr("y", function(d, i) {
            return (d.feature_type === 'CDS') ? 0 : lower_gene_chart_height/4;
        })
        .attr("height", function(d, i) {
            return (d.feature_type === 'CDS') ? lower_gene_chart_height : lower_gene_chart_height/2;
        });

    change_exon_plot(skip_utrs, exon_x_scale);

    // create strand-direction tooltip
    window.strand_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        return 'Strand Direction: ' + window.strand + '<br/>';
    });
    exon_track.call(window.strand_tip);

    // draw strand-direction arrow
    var a_s = window.strand == "-"? -1 : 1; //arrow direction
    var a_x = -5;  //arrow position on x-axis
    var a_y = lower_gene_chart_height/2.0; //arrow position on y-axis
    var a_w = 2; //arrow width
    var points = [[a_x+a_s*6, a_y], [a_x+a_s*1, a_y+a_w*3], [a_x+a_s*1, a_y+a_w], [a_x-a_s*9, a_y+a_w],
        [a_x-a_s*9, a_y-a_w], [a_x+a_s*1, a_y-a_w], [a_x+a_s*1, a_y-a_w*3]];
    exon_track.append("polygon")
            .attr("points", points.join(" "))
            .attr("fill", "steelblue")
            .attr("stroke", "black")
            .on('mouseover', window.strand_tip.show)
            .on('mouseout', window.strand_tip.hide);
}

function create_variants_plot(svg, exon_x_scale, coords) {
    var bounds = get_af_bounds(window.variants_for_graph);
    var min_af = bounds[0];
    var max_af = bounds[1];
    var variant_size_scale = d3.scale.log()
        .domain([min_af, max_af])
        .range([2, lower_gene_chart_height/3]);

    svg.selectAll('a.track_variant_link').data([]).exit().remove();
    svg.selectAll('a.track_variant_link').data(window.variants_for_graph).enter()
        .append("a")
        .attr('class', 'track_variant_link')
        .attr("xlink:href", function(d, i) { return "/variant/" + d.chrom + "-" + d.pos + "-" + d.ref + "-" + d.alt; })
        .attr("data-toggle", "tooltip")
        .attr('category', function(d) {
            return d.category;
        })
        .attr('variant_id', function(d) {
            return d.variant_id;
        })
        .on('mouseover', function(d) {
            $('#variant_' + d.variant_id).find('td').addClass('table_hover');
            window.variant_plot_tip.show(d);
        })
        .on('mouseout', function(d) {
            $('#variant_' + d.variant_id).find('td').removeClass('table_hover');
            window.variant_plot_tip.hide(d);
        })
        .append("ellipse")
        .attr("class", function(d) {
            return "track_variant " + d.category;
        })
        .style("opacity", 0.5)
        .style("visibility", function(d) {
            return (d[coords] == undefined) ? "hidden" : "visible";
        })
        .attr("cx",  function(d) {
            if (d[coords] == undefined) {
                // TODO: put these at the correct position, but hide them.
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        })
        .attr("cy", lower_gene_chart_height/2)
        .attr("rx", 2)
        .attr("ry", function(d, i) {
            if (!d.allele_freq) {
                return 0;
            } else {
                return variant_size_scale(d.allele_freq);
            }
        })
        // Workaround for exporting d3 to SVG (other delcaration in style.css).
        .attr('stroke', variant_colors)
        .attr('fill', variant_colors);
}

function add_variants_to_variants_plot() {
    var exon_track = d3.select('#gene_plot_container').select('#exon_track');
    var skip_utrs = ! $('#include_utrs_checkbox').is(':checked');
    var maximum_coding_position = get_maximum_coding_position(skip_utrs);
    var scale_type = $('.display_coverage_metric_buttons.active').attr('id').replace('display_coverage_', '').replace('_button', '');
    var chart_width = (scale_type === 'overview') ?
        gene_chart_width :
        Math.max(gene_chart_width, Math.min(maximum_coding_position, gene_chart_width*4));

    var exon_x_scale = d3.scale.linear()
        .domain([0, maximum_coding_position])
        .range([0, chart_width]);
    var coords = skip_utrs ? 'pos_coding_noutr' : 'pos_coding';

    create_variants_plot(exon_track, exon_x_scale, coords);
}

function variant_colors(d) {
    if (d.category == 'lof_variant') {
        return '#cd2932';
    } else if (d.category == 'missense_variant') {
        return '#a96500';
    } else if (d.category == 'synonymous_variant') {
        return '#157e28';
    }
}

function change_coverage_chart(coords, chart_width, exon_x_scale, metric) {
    var coords_start = coords.replace('pos_', 'start_');
    var coords_end = coords.replace('pos_', 'end_');

    var coverage_track = d3.select('#gene_plot_container').select('#coverage_svg')
        .attr("width", chart_width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .select('#coverage_track');

    var area = d3.svg.area()
        .x( function(d) { // Is .x() necessary?
            return exon_x_scale(d['pos']);
        }).y0( function(d) {
            return gene_chart_height;
        }).y1( function(d) {
            return (metric in d) ? y(d[metric]) : gene_chart_height;
        });

    //This is a hack, revealing that the surrounding code needs changes.
    var my_datum = [];
    _.each(window.coverage_stats, function(d) {
        if (d[coords_start] && d[coords_end]) {
            var o = [{}, {}, {}, {}];
            o[0]['pos'] = o[1]['pos'] = d[coords_start];
            o[2]['pos'] = o[3]['pos'] = 1 + d[coords_end];
            o[1][metric] = o[2][metric] = d[metric] || 0;
            o[0][metric] = o[3][metric] = 0;
            my_datum.push.apply(my_datum, o); //concatenate o to my_datum
        }
    });

    var max_cov = 1;
    if (metric == 'mean' || metric == 'median') {
        max_cov = d3.max(my_datum, function(d) { return d[metric]; });
    }
    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([gene_chart_height, 0]);

    var path = coverage_track.selectAll("path")
        .datum(my_datum)
        .transition()
        .duration(500)
        .attr("d", area)
        .style("fill", "steelblue")
        .attr('class', 'area');

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    coverage_track.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);
}

function change_exon_plot(skip_utrs, exon_x_scale) {
    var exon_track = d3.select('#exon_track');

    exon_track.selectAll("line.padded_exon")
        .transition()
        .duration(500)
        .style('visibility', function(d) {
            return (!skip_utrs || d.feature_type === 'CDS') ? 'visible' : 'hidden';
        })
        .attr("x1", function(d) {
            var x = exon_x_scale(get_coding_coordinate(d.start - EXON_PADDING, skip_utrs));
            return x;
        })
        .attr("x2", function(d) {
            var x = exon_x_scale(get_coding_coordinate(d.stop + EXON_PADDING, skip_utrs));
            return x;
        });

    $('#exon_track line.padded_exon').each(function (i, d) {
        if (d.getAttribute('x1') === '0' || d.getAttribute('x2') === '0') {
            d.style.visibility = 'hidden';
            $(d).hide();
            console.log(d);
        } else {
            //console.log(['nope!', d]);
        }
    });

    // plot exon rounded rects
    exon_track.selectAll("rect.exon_or_utr")
        .transition()
        .duration(500)
        .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(d.start, skip_utrs)); })
        .attr("width", function(d, i) {
            return (get_coding_coordinate(d.start, skip_utrs) === null) ? 0 : exon_x_scale(d.stop-d.start+1);
        })
        .style('visibility', function(d) {
            return (!skip_utrs || d.feature_type === 'CDS') ? 'visible' : 'hidden';
        });
}

function change_variant_plot(coords, exon_x_scale) {
    var exon_track = d3.select('#exon_track');
    exon_track.selectAll("a.track_variant_link")
        .transition()
        .duration(500)
        .selectAll('ellipse')
        .style("visibility", function(d) {
            return (d[coords] == undefined) ? "hidden" : "visible";
        })
        .attr("cx", function(d) {
            if (d[coords] == undefined) {
                // TODO: put these at the correct position, but hide them.
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        });
}

function change_plots(scale_type, metric, skip_utrs) {
    //scale_type is either "detail" or "overview"
    //metric is one of "mean", "median", "1", "5", "10", ... "50", "100"
    //skip_utrs is boolean

    var coords = skip_utrs ? 'pos_coding_noutr' : 'pos_coding';
    var maximum_coding_position = get_maximum_coding_position(skip_utrs);
    var chart_width = (scale_type == 'overview') ?
        gene_chart_width :
        Math.max(gene_chart_width, Math.min(gene_chart_width*4, maximum_coding_position*2));

    var exon_x_scale = d3.scale.linear()
        .domain([0, maximum_coding_position])
        .range([0, chart_width]);

    change_coverage_chart(coords, chart_width, exon_x_scale, metric);
    change_exon_plot(skip_utrs, exon_x_scale);
    change_variant_plot(coords, exon_x_scale);
}

function change_plots_with_values_from_page() {
    // Note: if this function runs immediately after a button is pressed, Bootstrap might not have toggled '.active' yet, so the results will be bad. Use setTimeout(..., 0)
    var skip_utrs = ! $('#include_utrs_checkbox').is(':checked');
    var scale_type = $('.display_coverage_metric_buttons.active').attr('id').replace('display_coverage_', '').replace('_button', '');
    var metric;
    if ($('.coverage_metric_buttons.active').attr('id') === 'covered_covmet_button') {
        metric = $('#over_x_select').val().replace('X', '');
    } else {
        metric = $("#average_select").val();
    }
    change_plots(scale_type, metric, skip_utrs);
}

function coverage_sum(key) {
    var total = 0;
    $.map(window.coverage_stats, function(entry) {
        total += entry[key];
    });
    return (total/window.coverage_stats.length).toPrecision(4);
}

$(document).ready(function() {
    build_the_graph();
    $("#coverage_plot_download").on('click', function() {
        download(set_plot_image('gene_plot_container', 0), $(this).attr('download'));
    });
    $("#exon_plot_download").on('click', function() {
        download(set_plot_image('gene_plot_container', 1), $(this).attr('download'));
    });
});

function download(data_uri, filename) {
    //This works in Chrome and Firefox and probably IE but not Safari.
    var link = document.createElement("a");
    link.download = filename;
    link.href = data_uri;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    delete link;
}

function set_plot_image(container, index) {
    //get svg element.
    var svg = $('#' + container).find('svg')[index];
    //get svg source.
    var serializer = new XMLSerializer();
    var source = serializer.serializeToString(svg);

    //add name spaces.
    if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)){
        source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
    }
    if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)){
        source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');
    }

    //add xml declaration
    source = '<?xml version="1.0" standalone="no"?>\r\n' + source;

    //convert svg source to URI data scheme.
    return "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(source);
}

function build_the_graph () {
    if ($(window).width() < 768) {
        $('#gene_plot_container').css('width', $(window).width() + "px");
    } else {
        $('#gene_plot_container').css('width', $(window).width()*10/12 + "px");
    }
    precalc_coding_coordinates_for_bins(window.coverage_stats); //Note: This modifies window.coverage_stats
    precalc_coding_coordinates(window.variants_for_graph); //Likewise for window.variants_for_graph

    // only show variants that have a coding coordinate
    window.variants_for_graph = _.filter(window.variants_for_graph, function(variant) {
        return variant.pos_coding != undefined;
    });

    // only show coding rects that have a coding coordinate
    window.coverage_stats = _.filter(window.coverage_stats, function(d) {
        return d.start !== undefined && d.end !== undefined;
    });
    $('#avg_coverage').html(coverage_sum('mean'));
    $('#avg_coverage_x').html(coverage_sum('30')*100 + '%');

    if (window.coverage_stats != null) {
        create_charts();
        if (window.variants_for_graph.length) {
            update_variants();
        }
        $('#loading_coverage').hide();
    } else {
        $('#gene_plot').hide();
        $('#not_covered').show();
    }

    // Listen for changes to the coverage plot
    $('.coverage_metric_buttons').change(function () {
        $('.coverage_subcat_selectors').hide();
        if ($(this).attr('id') === 'covered_covmet_button') {
            $('#over_x_select_container').show();
        } else {
            $('#average_select_container').show();
        }
        setTimeout(change_plots_with_values_from_page, 0);
    });
    $('#average_select').change(function () {
        $('#avg_coverage_type').html($(this).val());
        $('#avg_coverage').html(coverage_sum($(this).val()));
        setTimeout(change_plots_with_values_from_page, 0);
    });
    $('#over_x_select, #include_utrs_checkbox, .display_coverage_metric_buttons').change(function () {
        setTimeout(change_plots_with_values_from_page, 0);
    });
}
