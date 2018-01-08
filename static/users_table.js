$(function() {
    window.model = window.model || {};
    
    var search_timer = 0;
    var search_column = function(column, value) {
        return function() {
            if (column.search !== value) {
                column.search(value).draw();
            }
        }
    }

    var columns = [
            {title: 'User', name: 'username', data: 'username'}, 
            {title: 'E-mail', name: 'email', data: 'email'},
            {title: 'Enabled API', name: 'enabled_api', data: 'enabled_api'}
        ]

    $('#users_table').append('<tfoot><tr></tr></tfoot>');
    for (var i = 0; i < columns.length; i++) {
        $('#users_table tfoot tr').append('<th><input type="text" placeholder="Search ' + columns[i].title + '"/></th>')
    }

    window.model.tbl = $('#users_table').DataTable({
        serverSide: true,
        processing: true,
        deferRender: true,
        paging: true,
        pagingType: 'full',
        pageLength: 100,
        searching: true,
        scrollX: true,
        dom: 'ltipr',
        ajax: {
            url: window.model.url_prefix + 'administration/users',
            type: 'POST',
            data: function(args) {
                return { args: JSON.stringify(args) };
            } 
        },
        columns: columns
    });

    window.model.tbl.columns().every( function () {
        var that = this;
        $( 'input', this.footer() ).on( 'keyup change', function () {
            clearTimeout(search_timer);
            search_timer = setTimeout(search_column(that, this.value), 300);
        } );
    } );
});
