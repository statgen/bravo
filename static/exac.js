/*
 *
 *
 * Coding Coordinates
 *
 *
 */

var EXON_PADDING = 75;

/*
    The following methods are for working with "Coding Coordinates",
    a coordinate space that we use to plot data in the coding regions of a transcript.

    Conceptually: suppose you lined up all the coding regions in a transcript, with some padding on each side,
    then plotted any variants that overlap. The position of a variant on this line is the coding position -
    this will obviously differ from the actual genomic or cds coordinates.

    Random notes:
        - coding coordinates have no concept of a gene - they are solely a property of a transcript
        - should probably have a map between coding coodinates and protein position
 */
window.get_position_mapping = _.memoize(function(skip_utrs) {
    // Uses window.transcript.exons
    // Returns like [
    //    {real_start: 36649935, scaled_start: 0, length: 178},
    //    {real_start: 36650909, scaled_start: 179, length: 212}
    // ]

    var good_feature_types = skip_utrs ? ['CDS'] : ['CDS', 'UTR'];
    var exons = _.filter(window.transcript.exons, function(exon) {
        return _.contains(good_feature_types, exon.feature_type);
    });
    if (exons.length === 0) {
        exons = window.transcript.exons;
    }
    if (exons.length === 0) {
        return [];
    }

    // overlap b/w parts that touch is fine.
    exons = _.sortBy(exons, "start");
    var pos_mapping = [{
        real_start: exons[0].start - EXON_PADDING,
        scaled_start: 0,
        length: exons[0].stop - exons[0].start + EXON_PADDING*2
    }];
    for (var i=1; i<exons.length; i++) {
        var gap_between_exons = Math.min(exons[i].start - exons[i-1].stop, EXON_PADDING*2);
        var length_of_previous_exon = exons[i-1].stop - exons[i-1].start + 1; //Not sure about this +1.
        var scaled_start_of_previous_exon = pos_mapping[pos_mapping.length-1].scaled_start;
        pos_mapping.push({
            real_start: exons[i].start - EXON_PADDING,
            scaled_start: length_of_previous_exon + gap_between_exons + scaled_start_of_previous_exon,
            length: exons[i].stop - exons[i].start + EXON_PADDING*2
        });
    }
    console.log(pos_mapping);
    return pos_mapping;
});

window.get_coding_coordinates = function(positions, skip_utrs) {
    return positions.map(function(pos) {
        return window.get_coding_coordinate(pos, skip_utrs);
    })
};

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

window.get_coding_coordinate_params = function(skip_utrs) {
    //Seems to calculate the combined length of the exons including the padding between them.
    //With skip_utrs, each exon has EXON_PADDING on each side, meaning 2*EXON_PADDING between exons
    var ret = {};

    var pos_mapping = window.get_position_mapping(skip_utrs);
    ret.num_exons = pos_mapping.length;
    if (ret.num_exons === 0) {
        ret.size = 0;
    } else {
        ret.size = pos_mapping[pos_mapping.length-1].length + pos_mapping[pos_mapping.length-1].scaled_start - pos_mapping[0].scaled_start;
    }
    return ret;
};

window.precalc_coding_coordinates = function(objects, data_is_binned) {
    //Note: this modifies `objects`.
    var keys_to_map = data_is_binned ? ['start', 'end'] : ['pos'];
    keys_to_map.forEach(function(key) {
        orig_positions = _.map(objects, function(o) { return o[key] });
        new_positions = get_coding_coordinates(orig_positions, false);
        _.each(objects, function(o, i) {
            o[key + '_coding'] = new_positions[i];
        });
        new_positions = get_coding_coordinates(orig_positions, true);
        _.each(objects, function(o, i) {
            o[key + '_coding_noutr'] = new_positions[i];
        });
    });
};




/*
 *
 *
 * Other Stuff
 *
 *
 */

quality_chart_margin = {top: 10, right: 30, bottom: 45, left: 65};
quality_chart_height = 250 - quality_chart_margin.top - quality_chart_margin.bottom;
quality_chart_width = 300 - quality_chart_margin.left - quality_chart_margin.right;
xoffset = 40;
yoffset = 55;

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

function draw_region_coverage(raw_data, metric, ref) {
    // pjvh removed a bunch of functionality from this.  There's some useless code left behind.
    // If this function gets a single base, it draws the full distribution.
    // If it receives multiple bases, it draws a coverage graph, letting the user select mean, median, % > X
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

function update_variants() {
    var category_buttons = $('.consequence_display_buttons.active');
    if (category_buttons.length === 1) {
        var category = category_buttons.attr('id').replace('consequence_', '').replace('_variant_button', '');
        var filter = $('#filtered_checkbox').is(":checked") ? '[filter_status]' : '[filter_status="PASS"]';
        $('[category]').hide();
        if (category == 'other') {
            $('[category]' + filter).show();
        } else if (category == 'missense') {
            $('[category=missense_variant]' + filter).show();
        }
        $('[category=lof_variant]' + filter).show();
        if ($('tr[style!="display: none;"]').length == 1) {
            $('#variants_table_empty').show();
            $('#variants_table_container').hide();
        } else {
            $('#variants_table_empty').hide();
            $('#variants_table_container').show();
        }
    }
    $(document).trigger("just_updated_variants");
}

function get_af_bounds(data) {
    // Removing AC_Adj = 0 cases
    var min_af = d3.min(data, function(d) {
        if (d.allele_freq > 0) {
            return d.allele_freq;
        } else {
            return 1;
        }
    });
    // Should this be 1?
    var max_af = d3.max(data, function(d) { return d.allele_freq; });
    return [min_af, max_af];
}

total_width = $(window).width() < 768 ? $(window).width() : $(window).width() * 10 / 12;

gene_chart_margin = {top: 10, right: 10, bottom: 5, left: 30};
if ($(window).width() < 768) {
    gene_chart_margin.left = 10;
}
gene_chart_margin_lower = {top: 5, right: gene_chart_margin.right, bottom: 5, left: gene_chart_margin.left},
    gene_chart_width = total_width - gene_chart_margin.left - gene_chart_margin.right;

lower_gene_chart_height = 50 - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom,
    gene_chart_height = 300 - gene_chart_margin.top - gene_chart_margin.bottom - lower_gene_chart_height - gene_chart_margin_lower.top - gene_chart_margin_lower.bottom;

function change_track_chart_variant_size(variant_data, change_to, container) {
    var svg_outer = d3.select(container).select('#track');

    var variant_size_scale;
    var bounds = get_af_bounds(variant_data);
    var min_af = bounds[0];
    var max_af = bounds[1];
    if (change_to) {
        variant_size_scale = d3.scale.log()
            .domain([min_af, max_af])
            .range([lower_gene_chart_height / 3, 2]);
    } else {
        variant_size_scale = d3.scale.log()
            .domain([min_af, max_af])
            .range([2, lower_gene_chart_height / 3]);
    }
    svg_outer.selectAll("a")
        .selectAll("ellipse")
        .transition()
        .duration(500)
        .attr("ry", function(d, i) {
            if (!d.allele_freq) {
                return 0;
            } else {
                return variant_size_scale(d.allele_freq);
            }
        });
}

function memorySizeOf(obj) {
    var bytes = 0;

    function sizeOf(obj) {
        if(obj !== null && obj !== undefined) {
            switch(typeof obj) {
            case 'number':
                bytes += 8;
                break;
            case 'string':
                bytes += obj.length * 2;
                break;
            case 'boolean':
                bytes += 4;
                break;
            case 'object':
                var objClass = Object.prototype.toString.call(obj).slice(8, -1);
                if(objClass === 'Object' || objClass === 'Array') {
                    for(var key in obj) {
                        if(!obj.hasOwnProperty(key)) continue;
                        sizeOf(obj[key]);
                    }
                } else bytes += obj.toString().length * 2;
                break;
            }
        }
        return bytes;
    };

    function formatByteSize(bytes) {
        if(bytes < 1024) return bytes + " bytes";
        else if(bytes < 1048576) return(bytes / 1024).toFixed(3) + " KiB";
        else if(bytes < 1073741824) return(bytes / 1048576).toFixed(3) + " MiB";
        else return(bytes / 1073741824).toFixed(3) + " GiB";
    };

    return formatByteSize(sizeOf(obj));
}

// Adapted from http://jsfiddle.net/terryyounghk/KPEGU/
function exportTableToCSV($table, filename) {

    var $rows = $table.find('tr:has(td,th)[style!="display: none;"]'),

        // Temporary delimiter characters unlikely to be typed by keyboard
        // This is to avoid accidentally splitting the actual contents
        tmpColDelim = String.fromCharCode(11), // vertical tab character
        tmpRowDelim = String.fromCharCode(0), // null character

        // actual delimiter characters for CSV format
        colDelim = '","',
        rowDelim = '"\r\n"',

        // Grab text from table into CSV formatted string
        csv = '"' + $rows.map(function (i, row) {
            var $row = $(row),
                $cols = $row.find('td,th').not('.omit_csv');

            return $cols.map(function (j, col) {
                var $col = $(col),
                    text = $col.text();

                return text.replace('"', '""').replace(/\s+/g, " ").replace(/^\s+/, "").replace(/\s+$/, ""); // escape double quotes

            }).get().join(tmpColDelim);

        }).get().join(tmpRowDelim)
            .split(tmpRowDelim).join(rowDelim)
            .split(tmpColDelim).join(colDelim) + '"',

        // Data URI
        csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

    $(this)
        .attr({
        'download': filename,
            'href': csvData,
            'target': '_blank'
    });
}
function pad_2(number) { return (number < 10 ? '0' : '') + number; }

function date_format(date) {
     return date.getFullYear() + '_' +
         pad_2(date.getMonth()+1) + '_' +
         pad_2(date.getDate()) + '_' +
            pad_2(date.getHours()) + '_' +
            pad_2(date.getMinutes()) + '_' +
            pad_2(date.getSeconds()) ;
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