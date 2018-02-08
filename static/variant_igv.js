$(function() {
    window.model = window.model || {};

    var variant_id = window.variant.variant_id;
    var locus = window.variant.chrom + ":" + (window.variant.pos - 100) + "-" + (window.variant.pos + 100);

    var options = {
        locus: locus,
        showCommandBar: true,
        showKaryo: false,
        reference: {
           id: "hg38",
           fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa",
           indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.fai"
        }
    };

    var XHR = $.getJSON(window.model.url_prefix + 'variant/' + variant_id + '/reads');

    XHR.done( function(samples) {
        var n_hets = 0, n_homs = 0;
        for (var i = 0; i < samples.names.length; i++) {
            if (samples.names[i].startsWith("het-")) {
                n_hets += 1;
            } else if (samples.names[i].startsWith("hom-")) {
                n_homs += 1;
            }
        }

        var div = $("<div></div>").appendTo("#reads-het-samples :first-child");

        if (n_hets > 0) {
            div.addClass("btn-group").attr("role", "group");
        } else {
            div.addClass("alert alert-info").attr("role", "alert")
                .text("No heterozyhous individuals are present.");
        }

        div = $("<div></div>").appendTo("#reads-hom-samples :first-child");
        if (n_homs > 0) {
            div.addClass("btn-group").attr("role", "group");
        } else {
            div.addClass("alert alert-info").attr("role", "alert")
                .text("No homozygous individuals are present.");
        }

        var buttons = null;
        for (var i = 0; i < samples.names.length; i++) {
            if (samples.names[i].startsWith("het-")) {
                buttons = $("#reads-het-samples");
            } else if (samples.names[i].startsWith("hom-")) {
                buttons = $("#reads-hom-samples");
            } else {
                continue;
            }

            $("<button></button>").appendTo(buttons.find(".btn-group"))
                .addClass("btn btn-default")
                .attr("type", "button")
                .attr("id", samples.names[i])
                .click(function() {
                    if (!$(this).hasClass("active")) {
                        var track = {
                            id: this.id,
                            name: (this.id.startsWith("het-") ? "Heterozygous" : "Homozygous") + " Individual #" + this.id.split("-")[1],
                            indexed: true,
                            format: "bam",
                            type: "alignment",
                            colorBy: "strand",
                            url: window.variant.variant_id + "/" + this.id + ".bam"
                        };
                        igv_browser_instance.loadTrack(track);
                        $(this).addClass("active");
                        $(this).find(".glyphicon")
                            .removeClass("glyphicon-unchecked")
                            .addClass("glyphicon-check");
                    } else {
                        $(this).removeClass("active");
                        $(this).find(".glyphicon")
                            .removeClass("glyphicon-check")
                            .addClass("glyphicon-unchecked")
                        var trackPanelRemoved;
                        for (var j = 0; j < igv_browser_instance.trackViews.length; j++) {
                            if (this.id == igv_browser_instance.trackViews[j].track.id) {
                                trackPanelRemoved = igv_browser_instance.trackViews[j];
                                break;
                            }
                        }
                        if (trackPanelRemoved) {
                            igv_browser_instance.trackViews.splice(j, 1);
                            igv_browser_instance.trackContainerDiv.removeChild(trackPanelRemoved.trackDiv);
                            igv_browser_instance.fireEvent('trackremoved', [trackPanelRemoved.track]);
                        }
                    }
                })
                .append("<span class=\"glyphicon glyphicon-unchecked\"></span> Individual #" + samples.names[i].split("-")[1]);
        }

        $("#reads-spinner").addClass("hidden");
        $("#reads-het-samples-header").removeClass("hidden");
        $("#reads-het-samples").removeClass("hidden");
        $("#reads-hom-samples-header").removeClass("hidden");
        $("#reads-hom-samples").removeClass("hidden");
        $("#reads-igv").removeClass("hidden");

        var igv_browser_instance = igv.createBrowser($("#reads-igv-instance"), options);
        $(".igv-ideogram-content-div").hide();

        igv_browser_instance.on("trackremoved", function(track) {
            var button = $("button#" + track.id);
            button.removeClass("active");
            button.find('.glyphicon').removeClass('glyphicon-check').addClass("glyphicon-unchecked");
        });

        if (n_hets > 0) {
            $("[id^='het-']")[0].click();
        } else if (n_homs > 0) {
            $("[id^='hom-']")[0].click();
        }
    });
});
