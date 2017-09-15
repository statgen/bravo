create_summary_table();
create_coverage_plot();
$(function() {
    transcripts_plot.create();
    create_variant_plot();
    create_variant_table();
    create_pos_plot();
});
