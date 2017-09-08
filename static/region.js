
create_summary_table();
create_coverage_plot();
$(function() {
    populate_variant_table_filters();
    create_variant_plot();
    create_variant_table();
    initiate_mouse_guide();
});
