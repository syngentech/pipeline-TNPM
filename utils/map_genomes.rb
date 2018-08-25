#!/usr/bin/env ruby
require 'yaml'
dict = YAML.load_file(File.expand_path('map_genomes.dict.yaml', __dir__))

STDIN.each do |line|
  if line[0..1] != '##'
    col = line.chomp.split("\t")
    next unless dict[col[0]]

    col[0] = dict[col[0]]
    line = col.join("\t")
  end

  puts line
end
