use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

=begin listdb

$registry->load_registry_from_db(
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous'
);

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

foreach my $db_adaptor (@db_adaptors) {
   my $db_connection = $db_adaptor->dbc();

   printf(
      "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
      $db_adaptor->species(), $db_adaptor->group(),
      $db_connection->dbname(), $db_connection->host(),
      $db_connection->port()
   );
}

=end listdb

=cut

my $dbname = $ARGV[0];

$registry->load_registry_from_db(
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
   -dbname => $dbname
); 

my $genes_adaptor = $registry->get_adaptor('human', 'core', 'gene');
my @genes = @{ $genes_adaptor->fetch_all };

foreach my $gene (@genes) {
   printf("%s %s\n", $gene->stable_id, $gene->canonical_transcript->stable_id)
}
