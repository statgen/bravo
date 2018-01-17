$(function() {
    window.model = window.model || {};

    var options = {
            locus: "22:16,389,347-16,389,547",
            showCommandBar: true,
            showKaryo: false,

            reference: {
               id: "hg38"
            },

            tracks: [
               {
                   name: "test_track",
                   url: "22-16389447-A-G/test.bam"
               }
            ]
        };

    var div = $("#igvDiv")[0]

    var browser = igv.createBrowser(div, options);    
});
