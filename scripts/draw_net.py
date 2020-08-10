#!/usr/bin/python3

import argparse
import os
import re
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Draw Net')
parser.add_argument('input_file_names', nargs='+')
parser.add_argument('-d', '--m1_direction', default='horizontal')
args = parser.parse_args()


class Box:
    def __init__(self, layer, lx, hx, ly, hy):
        self.layer = int(layer)
        self.lx = int(lx)
        self.hx = int(hx)
        self.ly = int(ly)
        self.hy = int(hy)

    def __repr__(self):
        return 'box(l={}, x=({}, {}), y=({}, {}))'.format(self.layer, self.lx, self.hx, self.ly, self.hy)

    def w(self):
        return self.hx - self.lx

    def h(self):
        return self.hy - self.ly

def intersect(b1, b2):
    lx = max(b1.lx, b2.lx)
    ly = max(b1.ly, b2.ly)
    hx = min(b1.hx, b2.hx)
    hy = min(b1.hy, b2.hy)
    return lx<=hx and ly<=hy


for file_name in args.input_file_names:
    pinNames = []
    pinAccessBoxes = []
    oriRGs = []
    paths = []
    pickedRGs = []
    px = []
    py = []

    reHeader = re.compile('Net (.*) \(idx = (\d+)\) with (\d+) pins')
    rePin = re.compile('pin (\d+) p.*')
    reGuide = re.compile('(\d+) route guides')
    reBox = re.compile('box\(l=(\d+), x=\((-?\d+), (-?\d+)\), y=\((-?\d+), (-?\d+)\)\)')
    with open(file_name) as file:
        result = None
        while not result:
            line = file.readline()
            result = reHeader.search(line)
            netName, netIdx, numPins = result.group(1), int(result.group(2)), int(result.group(3))

        # pin
        for _ in range(numPins):
            result = rePin.search(next(file))
            pinNames.append(result.group(1))
            result = reBox.search(next(file))
            pinAccessBoxes.append(Box(result.group(1), result.group(2), result.group(3), result.group(4), result.group(5)))

        # original guide
        numGuides = int(reGuide.search(next(file)).group(1))
        for _ in range(numGuides):
            result = reBox.search(file.readline())
            oriRGs.append(Box(result.group(1), result.group(2), result.group(3), result.group(4), result.group(5)))

        # paths
        numPaths = int(next(file).split(' ')[0])
        for _ in range(numPaths):
            result = reBox.search(file.readline())
            paths.append(Box(result.group(1), result.group(2), result.group(3), result.group(4), result.group(5)))
        
        # picked guide
        numPick = int(next(file).split(' ')[0])
        for _ in range(numPick):
            result = reBox.search(file.readline())
            pickedRGs.append(Box(result.group(1), result.group(2), result.group(3), result.group(4), result.group(5)))

        # pnet lines
        next(file)
        line = next(file)
        while len(line) > 2:
            result = reBox.search(line)
            box = Box(result.group(1), result.group(2), result.group(3), result.group(4), result.group(5))
            px.append([box.lx])
            py.append([box.ly])
            line = next(file)
            result = reBox.search(line)
            box = Box(result.group(1), result.group(2), result.group(3), result.group(4), result.group(5))
            px[-1].append(box.lx)
            py[-1].append(box.ly)
            line = next(file)

    # num of layers
    numLayers = 0
    bbox = Box(0, 100000000, 0, 100000000, 0)
    layerSet = set()
    for box in pinAccessBoxes:
        numLayers = max(box.layer, numLayers)
        layerSet.add(box.layer)
    for box in paths:
        numLayers = max(box.layer, numLayers)
        layerSet.add(box.layer)
    for box in pickedRGs:
        numLayers = max(box.layer, numLayers)
        bbox.lx = min(box.lx, bbox.lx)
        bbox.ly = min(box.ly, bbox.ly)
        bbox.hx = max(box.hx, bbox.hx)
        bbox.hy = max(box.hy, bbox.hy)
        layerSet.add(box.layer)
    numLayers += 1
    print('# layers is {}'.format(numLayers))

    # define box plotting
    cmap = plt.get_cmap('plasma')
    colors = [cmap(i / (numLayers - 1)) for i in range(0, numLayers)]
    hatches = ['--', '||', '//', '\\', '++', 'xx', 'oo', 'OO', '..', '**']
    horiOffset = 0 if args.m1_direction == 'horizontal' else 1
    labels = ['M{} '.format(i+1) + ('horizontal' if i % 2 == horiOffset else "vertical") for i in range(0, numLayers)]
    def plotBox(box):
        plt.gca().add_patch(plt.Rectangle((box.lx, box.ly), box.w(), box.h(), lw=1,
                                        edgecolor='k', facecolor=colors[box.layer], alpha=0.3, hatch=hatches[box.layer], label=labels[box.layer]))
    def plotBoxNoLabel(box):
        plt.gca().add_patch(plt.Rectangle((box.lx, box.ly), box.w(), box.h(), lw=1,
                                        edgecolor='k', facecolor=colors[box.layer], alpha=0.3, hatch=hatches[box.layer]))
    def plotPaths(box, label):
        if not label:
            plt.gca().add_patch(plt.Rectangle((box.lx, box.ly), box.w(), box.h(), lw=2,
                                        edgecolor='black', facecolor=colors[box.layer], alpha=0.6, hatch=hatches[box.layer], label=labels[box.layer]))
        else:
            plt.gca().add_patch(plt.Rectangle((box.lx, box.ly), box.w(), box.h(), lw=2,
                                        edgecolor='black', facecolor=colors[box.layer], alpha=0.6, hatch=hatches[box.layer]))


    # plot boxes (in the order of layers to make legend ordered)
    plt.figure(figsize=(15, 30))
    # plt.axes()
    plt.subplot(122)
    layerl = [0,0,0,0,0,0,0,0,0]
    for box in pinAccessBoxes:
        box.lx -= 20
        box.hx += 20
        box.ly -= 20
        box.hy += 20
        plotBoxNoLabel(box)
    for box in pickedRGs:
        plotBoxNoLabel(box)
    for box in paths:
        box.lx -= 100
        box.hx += 100
        box.ly -= 100
        box.hy += 100
        plotPaths(box, layerl[box.layer])
        layerl[box.layer] = 1
    plt.legend()

    # anno pins
    for pinIdx in range(numPins):
        pinX = pinAccessBoxes[pinIdx].lx
        pinY = pinAccessBoxes[pinIdx].ly
        anno = 'pin {} M{}'.format(pinNames[pinIdx], pinAccessBoxes[pinIdx].layer+1)
        plt.text(pinX, pinY, anno, color='k')
    
    # pnet lines
    for i in range(len(px)):
        plt.plot(px[i],py[i],c='r')

    # format
    plt.xlabel('X (DBU)')
    plt.ylabel('Y (DBU)')
    plt.title('{} with {} pins and {} paths'.format(netName, numPins, numPaths))
    plt.axis('square')
    # plt.axis('scaled')

    # another
    plt.subplot(121)
    layerl = [0,0,0,0,0,0,0,0,0]
    again = True
    for box in oriRGs:
        if intersect(box, bbox) and (box.layer in layerSet):
            if (layerl[box.layer] > 0):
                plotBoxNoLabel(box)
            else:
                plotBox(box)
                layerl[box.layer] = 1
            again = False
    # if again:
    #     numLayers = 9
    #     labels = ['M{} '.format(i+1) + ('horizontal' if i % 2 == horiOffset else "vertical") for i in range(0, numLayers)]
    #     colors = [cmap(i / (numLayers - 1)) for i in range(0, numLayers)]
    #     for box in oriRGs:
    #         if layerl[box.layer] > 0:
    #             plotBoxNoLabel(box)
    #         else:
    #             plotBox(box)
    #             layerl[box.layer] = 1
    plt.legend()
    plt.axis('square')
    # plt.axis('scaled')

    # save
    base_name, ext_name = os.path.splitext(os.path.basename(file_name))
    # plt.savefig('{}.pdf'.format(base_name), bbox_inches='tight')
    # plt.savefig('{}.png'.format(base_name), bbox_inches='tight')
    plt.show()