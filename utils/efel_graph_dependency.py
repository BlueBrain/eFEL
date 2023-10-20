'''graph the dependency graph'''

import argparse
import logging
import os

from collections import namedtuple, defaultdict

Feature = namedtuple('Feature', 'library name dependencies')
Dependency = namedtuple('Dependency', 'feature wildcard')


def create_feature(line):
    '''Convert a line from the Dependency.txt to a feature

    LibV2:AP_duration_half_width_change	#LibV2:AP_duration_half_width
    LibV2:E6	#LibV1:AP_amplitude;APWaveForm
    '''
    tokenized_line = line.strip().split()

    library, name = '', tokenized_line[0]
    if ':' in name:
        library, name = name.strip().split(':')

    dependencies = []
    if len(tokenized_line) > 1:
        dependencies.extend(create_dependency(d) for d in tokenized_line[1:])

    return Feature(library=library, name=name, dependencies=dependencies)


def create_dependency(line):
    '''parse a dependency, and create a Dependency'''
    line = line.strip()
    assert line[0] == '#', '"%s": Not a feature dependency, must start w/ #' % line
    feature, wildcard = line[1:], ''
    if ';' in feature:
        feature, wildcard = feature.split(';')
    return Dependency(feature=feature, wildcard=wildcard)


def get_feature_library(dependency_file):
    '''parse the `dependency_file` and gather all the features and their dependencies'''
    feature_library = defaultdict(dict)
    with open(dependency_file) as fd:
        for line in fd:
            feature = create_feature(line)
            feature_library[feature.library][feature.name] = feature
            #  defaults are the last to be loaded
            feature_library['Default'][feature.name] = feature

    def find_feature(dependency):
        '''find the feature in the current feature_library'''
        if ':' in dependency.feature:
            library, name = dependency.feature.split(':')
        else:
            library, name = 'Default', dependency.feature
        feature = feature_library[library][name]
        return Dependency(feature=feature, wildcard=dependency.wildcard)

    #  resolve features pointers
    for library in feature_library.keys():
        if 'Default' == library:
            continue
        for feature in feature_library[library].values():
            for i in range(len(feature.dependencies)):
                feature.dependencies[i] = find_feature(feature.dependencies[i])

    return feature_library


COLORS = ('blue', 'crimson', 'forestgreen', 'dodgerblue3', 'aquamarine4',
          )


def output_graphivz(feature_library, filename, draw_dependencies=False,
                    output_dot=False):
    '''create a PNG at filename of the dependencies'''
    try:
        import pygraphviz
    except ImportError:
        import sys
        sys.exit('Need to have the "pygraphviz" package installed')
    G = pygraphviz.AGraph(name='Dependencies', directed=True)
    # below attributes are for the aesthetics
    G.graph_attr['overlap'] = 'false'
    G.graph_attr['rankdir'] = 'LR'
    G.graph_attr['ranksep'] = '2.5'
    G.graph_attr['nodesep'] = '0'


    for library, color in zip(sorted(feature_library.keys()), COLORS):
        if 'Default' == library:
            continue
        G.add_node(library, shape='circle', color=color)
        sg = G.add_subgraph(name=library, style="")
        for name, feature in feature_library[library].items():
            full_name = library + ':' + name
            sg.add_node(full_name, label=name, shape='box', color=color)
            sg.add_edge(library, full_name, color=color)
            if draw_dependencies:
                for dep in feature.dependencies:
                    target = dep.feature.library + ':' + dep.feature.name
                    G.add_edge(full_name, target, color='yellowgreen')

    if output_dot:
        G.write(os.path.splitext(filename)[0] + '.dot')
    G.draw(filename, prog='dot')


def get_parser():
    '''return the argument parser'''
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True,
                        help='Input Dependency.txt file')
    parser.add_argument('-g', '--graph',
                        help='Output graph png file location')
    parser.add_argument('--graph-deps', action='store_true',
                        help='If graphing, draw the links to the dependencies')
    parser.add_argument('--dot', action='store_true',
                        help='If graphing, output the dot file')
    parser.add_argument('-d', '--dump', action='store_true',
                        help='Dump the resolved dependency tree')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='-v for INFO, -vv for DEBUG')
    return parser


def main(args):
    '''main'''
    logging.basicConfig(level=(
        logging.WARNING, logging.INFO, logging.DEBUG)[min(args.verbose, 2)])

    feature_library = get_feature_library(args.input)
    if args.dump:
        from pprint import pprint
        pprint(dict(feature_library))
    if args.graph:
        output_graphivz(feature_library, args.graph, args.graph_deps, args.dot)


if __name__ == '__main__':
    PARSER = get_parser()
    main(PARSER.parse_args())
