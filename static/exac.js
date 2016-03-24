function update_variants() {
    var category_buttons = $('.consequence_display_buttons.active');
    if (category_buttons.length === 1) {
        var category = category_buttons.attr('id').replace('consequence_', '').replace('_variant_button', '');
        $('[category]').hide();
        if (category == 'other') {
            $('[category]').show();
        } else if (category == 'missense') {
            $('[category=missense_variant]').show();
        }
        $('[category=lof_variant]').show();
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

function with_waiting_notice(f) {
    $('body').css('cursor', 'progress');
    $("#wait-modal").modal({"backdrop": "static"});
    setTimeout(function() {
        f();
        $('body').css('cursor', 'default');
        $("#wait-modal").modal('hide');
    }, 10);
}
