local class = require 'ext.class'
local template = require 'template'

local CLKernel = class()

function CLKernel:init(args)
	self.env = assert(args.env)
	self.name = args.name or 'kernel_'..tostring(self):sub(10)
	self.argsOut = args.argsOut
	self.argsIn = args.argsIn
	self.argBuffers = table()
		:append(self.argsOut)
		:append(self.argsIn)
		:map(function(arg) return arg.buf end)
	self.code = table{
		args.header or '',
		template([[
kernel void <?=self.name?>(
<?
local sep = ''
for _,arg in ipairs(self.argsOut or {}) do 
?>	<?=sep?>global <?=arg.type?>* <?=arg.name?>
<?
sep = ', '
end
for _,arg in ipairs(self.argsIn or {}) do
?>	<?=sep?>global const <?=arg.type?>* <?=arg.name?>
<?
sep = ', '
end
?>
) {
INIT_KERNEL();
<?=args.body?>
}
]], {self=self, args=args})
	}:concat'\n'
end

function CLKernel:__call(...)
	-- if we get a call request when we have no kernel/program, make sure to get one 
	if not self.kernel then
		self.program = require 'cl.program'{context=self.env.ctx, devices={self.env.device}, code=self.code}
		self.kernel = self.program:kernel(self.name, self.argBuffers:unpack())
	end

	self.env:clcall(self.kernel, ...)
end

return CLKernel
