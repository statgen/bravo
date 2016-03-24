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

function with_waiting_notice(f) {
    $('body').css('cursor', 'progress');
    $("#wait-modal").modal({"backdrop": "static"});
    setTimeout(function() {
        f();
        $('body').css('cursor', 'default');
        $("#wait-modal").modal('hide');
    }, 10);
}
