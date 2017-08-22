window.model = window.model || {};
window._debug = window._debug || {};

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

(function() {
    var autocomplete_bloodhound = new Bloodhound({
        datumTokenizer: Bloodhound.tokenizers.obj.whitespace('value'),
        queryTokenizer: Bloodhound.tokenizers.whitespace,
        identify: function(sugg) { return sugg.value; }, // maybe allows Bloodhound to `.get()`  objects
        remote: {
            url: '/api/autocomplete?query=%QUERY',
            wildcard: '%QUERY',
            rateLimitBy: 'throttle',
            rateLimitWait: 500,
            transform: function(data) {
                // Probably this function reveals that I don't understand Bloodhound.
                // But, I want my previous results to stay around while I keep typing.
                // If the string that's currently in the searchbox matches some string that's been suggested before, I want to see it!
                // This especially happens while I'm typing a chrom-pos-ref-alt.  If what I'm typing agrees with something being suggested, it shouldn't disappear!
                // So, I'm just adding everything to the local index. (Note: NOT localstorage.)
                // Bloodhound appears to perform deduping.
                autocomplete_bloodhound.add(data);
                return data;
            },
        },
        sorter: function(a, b) { return (a.value > b.value) ? 1 : -1; },
    });
    autocomplete_bloodhound.initialize();

    $(function() {
        $('.typeahead').typeahead({
            hint: false,
            highlight: true,
            minLength: 1,
        }, {
            name: 'autocomplete',
            source: autocomplete_bloodhound,
            display: 'value',
            limit: 10,
            templates: {
                suggestion: _.template("<div><%= value %></div>"),
                // Currently, variants and rsids do not give autocomplete results.
                // So, if `empty` is enabled, then they always say "No matches found", even though there ARE matches.
                // So, I'm removing this.  A bad solution, but better than nothing.
                //empty: "<div class='tt-empty-message'>No matches found.</div>"
            }
        });

        $('.typeahead').bind('typeahead:select', function(ev, suggestion) {
            window.location.href = '/awesome?query=' + suggestion.value;
        });
    });
})();


// convenience functions
function fmt(format) {
    var args = Array.prototype.slice.call(arguments, 1);
    return format.replace(/{(\d+)}/g, function(match, number) {
        return (typeof args[number] != 'undefined') ? args[number] : match;
    });
}

function group_thousands(x) {
    // group_thousands_space(1234567.7654321) -> "1 234 567.765 432 1"
    try {
        var thinspace = '\u202f';
        var parts = (''+x).split('.');
        var b = parts[0]; // before the decimal
        var neg = b[0] === '-';
        if (neg) b = b.substr(1);
        var L = b.length;
        var i = b.length;
        var bf = ''; // before the decimal, formatted
        while (i--) bf = (i===0?'':((L-i)%3?'':thinspace)) + b.charAt(i) + bf;
        if (neg) bf = '-' + bf;
        if (parts.length === 1) return bf;
        var a = parts[1]; // after the decimal
        L = a.length;
        i = 0;
        var af = ''; // after the decimal, formatted
        while (i < L) { af += (i===0?'':(i%3?'':thinspace)) + a.charAt(i); i++; }
        return bf + '.' + af;
    } catch(e) {
        console.error(['exception occurrend in group_thousands', e, x]);
        return x;
    }
}
function perc_sigfigs(d, n_sigfigs, n_left_of_decimal) {
    /* like proportion_sigfigs but outputs a percentage.
       left-pads with &nbsp; to line up `100%` with `0%`
       a leading 9 will never be rounded up, so the number of digits can never change. */
    try {
        if (typeof n_sigfigs === 'undefined') n_sigfigs = 2;
        if (typeof n_left_of_decimal === 'undefined') n_left_of_decimal = 3;
        var nobreakspace = '\u2007'; // '\u00A0';

        d = d*100;
        var parts = d.toString().replace(/^0+/, '').split('.');
        var decimal_offset = parts[0].length;
        if (parts[0].length > n_sigfigs) n_sigfigs = parts[0].length; /* I hate zero used to imply "dunno" */
        var onestring = parts.join('');
        onestring = sigfig(onestring, n_sigfigs);
        if (onestring.length < decimal_offset) { throw 'wat ' + onestring + ' ' + decimal_offset }
        var before_decimal = onestring.substr(0, decimal_offset);
        var after_decimal = onestring.substr(decimal_offset);
        if (before_decimal.length === 0) before_decimal = '0';
        before_decimal = leftpad(before_decimal, n_left_of_decimal, nobreakspace);
        var full = before_decimal + (after_decimal.length ? '.' + after_decimal : '');
        full = group_thousands(full);
        return full + '%'; // + ' (' + d.toString() + ')';
    } catch(e) {
        console.error(['exception occurrend in perc_sigfigs', e, [d, n_sigfigs, n_left_of_decimal]]);
        return d;
    }
}
function sigfig(numstring, k) {
    /* sigfigs("999499", 1) -> "9995" */
    var not_sig = ''
    if (numstring[0] === '9') not_sig = numstring.match(/^9*/)[0];
    if (numstring[0] === '0') not_sig = numstring.match(/^0*/)[0];
    var sig = numstring.substr(not_sig.length);
    var firstk = sig.substr(0, k);
    if (sig.length > k) { /* round */
        var rounder = sig[k];
        if (parseInt(rounder) >= 5) firstk = (parseInt(firstk) + 1).toString();
    }
    return not_sig + firstk;
}
function leftpad(str, len, padding) {
    str = String(str); len = len>>0; if (str.length >= len) return str; return padding.repeat(len - str.length) + str;
}
function rightpad(str, len, padding) {
    str = String(str); len = len>>0; if (str.length >= len) return str; return str + padding.repeat(len - str.length);
}
function fmt_annotation(anno) {
    return anno.replace('_variant', '').replace(/_/g, ' ').replace('utr', 'UTR').replace('3 prime', "3'").replace('5 prime', "5'").replace('nc ', "non-coding ");
}
