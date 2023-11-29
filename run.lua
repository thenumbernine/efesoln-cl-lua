#!/usr/bin/env luajit

-- global
cmdline = require 'ext.cmdline'(...)

local App = require 'app'
local app = App()
return app:run()
