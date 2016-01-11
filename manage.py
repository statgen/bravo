#!/usr/bin/env python2

from flask.ext.script import Manager
import exac

manager = Manager(exac.app)


@manager.command
def hello():
    print "hello"


@manager.command
def load_db():
    exac.load_db()


@manager.command
def load_base_coverage():
    exac.load_base_coverage()


@manager.command
def load_variants_file():
    exac.load_variants_file()


@manager.command
def load_gene_models():
    exac.load_gene_models()


@manager.command
def load_dbsnp_file():
    exac.load_dbsnp_file()


@manager.command
def create_cache():
    exac.create_cache()


@manager.command
def precalculate_metrics():
    exac.precalculate_metrics()


@manager.command
def precalculate_variant_consqequence_category():
    exac.precalculate_variant_consqequence_category()

@manager.command
def create_users():
    exac.create_users()

if __name__ == "__main__":
    manager.run()

