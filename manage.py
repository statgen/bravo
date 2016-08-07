#!/usr/bin/env python2

from flask_script import Manager
import exac

manager = Manager(exac.app)


@manager.command
def hello():
    print "hello"


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
def precalculate_metrics():
    exac.precalculate_metrics()


@manager.command
def precalculate_whether_variant_is_ever_missense_or_lof():
    exac.precalculate_whether_variant_is_ever_missense_or_lof()

@manager.command
def create_users():
    exac.create_users()

if __name__ == "__main__":
    manager.run()

