
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
    $('.frequency_display_buttons').change(function() {
        $('.frequency_displays').hide();
        var v = $(this).attr('id').replace('_button', '');
        $('#' + v + '_container').show();
    });

    $('#frequency_table').tablesorter({
        stringTo: 'bottom',
        sortList: [[4,1], [0,0]],
        headers: {
            4: {
                sorter: 'digit'
            }
        }
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
