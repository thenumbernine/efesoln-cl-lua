#!/usr/bin/env luajit
--[[
this should be a stand-alone tool
--]]
require 'ext'
local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local sdl = require 'ffi.sdl'
local ig = require 'ffi.imgui'
local ImGuiApp = require 'imguiapp'
local Mouse = require 'gui.mouse'
local quat = require 'vec.quat'
local vec3 = require 'vec.vec3'
local vec4 = require 'vec.vec4'
local vec4d = require 'vec-ffi.vec4d'
local glreport = require 'gl.report'
local glcall = require 'gl.call'
local GradientTex = require 'gl.gradienttex'
local Tex2D = require 'gl.tex2d'
local Tex3D = require 'gl.tex3d'
local GLProgram = require 'gl.program'
local EFESolver = require 'efe'

local mouse = Mouse()
local viewAngle = quat()
local viewDist = 2
	
local hsvTex
local volumeShader

local clipEnabled = ffi.new('bool[1]', true)
local rotateClip = ffi.new('int[1]', 0)

local clipInfos = range(4):map(function(i)
	local plane = vec4(0,0,0,0)
	plane[math.min(i,3)] = -1
	return {
		enabled = ffi.new('bool[1]', i==3),
		plane = plane,
	}
end)

local alpha = ffi.new('float[1]', 1.5e-1)
local alphaGamma = ffi.new('float[1]', 1)
local showGradTrace = ffi.new('bool[1]', false)
local showCurlTrace = ffi.new('bool[1]', false)


local App = class(ImGuiApp)

App.title = 'Einstein Field Equation Solver'

function App:initGL()
	App.super.initGL(self)

	self.solver = EFESolver{
		app = self,
		config = require 'config',
	}

	local hsvWidth = 256
	hsvTex = GradientTex(hsvWidth,
--[[ rainbow or heatmap or whatever
		{
			{0,0,0,0},
			{1,0,0,1/6},
			{1,1,0,2/6},
			{0,1,1,3/6},
			{0,0,1,4/6},
			{1,0,1,5/6},
			{0,0,0,6/6},
		},
--]]
-- [[ sunset pic from https://blog.graphiq.com/finding-the-right-color-palettes-for-data-visualizations-fcd4e707a283#.inyxk2q43
		table{
			vec3(22,31,86),
			vec3(34,54,152),
			vec3(87,49,108),
			vec3(156,48,72),
			vec3(220,60,57),
			vec3(254,96,50),
			vec3(255,188,46),
			vec3(255,255,55),
		}:map(function(c,i)
			return table(c/255):append{1}
		end),
--]]
		false)
	-- change to 2D so imgui can use it
	local data = ffi.new('unsigned char[?]', hsvWidth*4)
	gl.glGetTexImage(gl.GL_TEXTURE_1D, 0, gl.GL_RGBA, gl.GL_UNSIGNED_BYTE, data)
	hsvTex:unbind()
	hsvTex:delete()
	hsvTex = Tex2D{
		internalFormat = gl.GL_RGBA,
		width = hsvWidth,
		height = 1,
		format = gl.GL_RGBA,
		type = gl.GL_UNSIGNED_BYTE,
		data = data,
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		magFilter = gl.GL_LINEAR,
		wrap = {
			s = gl.GL_CLAMP_TO_EDGE,
			t = gl.GL_REPEAT,
		},
		generateMipmap = true,
	}

	volumeShader = GLProgram{
		vertexCode = [[
varying vec3 pos;
void main() {
	pos = gl_Vertex.xyz;
	gl_Position = ftransform();
}
]],
		fragmentCode = [[
varying vec3 pos;
uniform sampler3D volTex;
uniform sampler2D hsvTex;
uniform vec3 normal;
uniform float alpha;
uniform float alphaGamma;
uniform bool clipEnabled[4];
void main() {

	vec4 worldPos = gl_ModelViewMatrix * vec4(pos,1.);
	for (int i = 0; i < 4; ++i) {
		if (clipEnabled[i] && dot(worldPos, gl_ClipPlane[i]) < 0.) discard;
	}

	float value = texture3D(volTex, pos).r;
	vec4 voxelColor = vec4(texture2D(hsvTex, vec2(value, .5)).rgb, pow(alpha, alphaGamma));
	
	//calculate normal in screen coordinates
	vec4 n = gl_ModelViewProjectionMatrix * vec4(normal, 0.);
	//determine length of line through slice at its angle
	voxelColor.a /= -n.w;
	
	gl_FragColor = vec4(voxelColor.rgb, voxelColor.a * alpha);
}
]],
		uniforms = {
			volTex = 0,
			hsvTex = 1,
		},
	}

	glreport'here'

	gl.glEnable(gl.GL_DEPTH_TEST)
end

local leftShiftDown
local rightShiftDown 
function App:event(event, eventPtr)
	App.super.event(self, event, eventPtr)
	local canHandleMouse = not ig.igGetIO()[0].WantCaptureMouse
	local canHandleKeyboard = not ig.igGetIO()[0].WantCaptureKeyboard
	
	if event.type == sdl.SDL_MOUSEBUTTONDOWN then
		if event.button.button == sdl.SDL_BUTTON_WHEELUP then
			orbitTargetDistance = orbitTargetDistance * orbitZoomFactor
		elseif event.button.button == sdl.SDL_BUTTON_WHEELDOWN then
			orbitTargetDistance = orbitTargetDistance / orbitZoomFactor
		end
	elseif event.type == sdl.SDL_KEYDOWN or event.type == sdl.SDL_KEYUP then
		if event.key.keysym.sym == sdl.SDLK_LSHIFT then
			leftShiftDown = event.type == sdl.SDL_KEYDOWN
		elseif event.key.keysym.sym == sdl.SDLK_RSHIFT then
			rightShiftDown = event.type == sdl.SDL_KEYDOWN
		elseif canHandleKeyboard and event.type == sdl.SDL_KEYDOWN then
			if event.key.keysym.sym == sdl.SDLK_UP then
				self.solver.displayVarPtr[0] = math.max(0, self.solver.displayVarPtr[0] - 1)
				self.solver:refreshDisplayVarKernel()
			elseif event.key.keysym.sym == sdl.SDLK_DOWN then
				self.solver.displayVarPtr[0] = math.min(#self.solver.displayVars-1, self.solver.displayVarPtr[0] + 1)
				self.solver:refreshDisplayVarKernel()
			elseif event.key.keysym.sym == sdl.SDLK_SPACE then
				self.updateMethod = not self.updateMethod
			elseif event.key.keysym.sym == ('u'):byte() then
				self.updateMethod = 'step'
			elseif event.key.keysym.sym == ('r'):byte() then
				print'resetting...'
				self.solver:resetState()
				self.updateMethod = nil
			end
		end
	end
end

App.updateMethod = nil

function App:update()
	if self.updateMethod then
		if self.updateMethod == 'step' then 
			print('performing single step...')
			self.updateMethod = nil 
		end
		
		self.solver:update()
	end

	local canHandleMouse = not ig.igGetIO()[0].WantCaptureMouse
	if canHandleMouse then 
		mouse:update()
	end
	if mouse.leftDragging then
		if leftShiftDown or rightShiftDown then
			if rotateClip[0] == 0 then
				viewDist = viewDist * math.exp(10 * mouse.deltaPos[2])
			else
				local clipPlane = clipInfos[rotateClip[0]].plane
				clipPlane[4] = clipPlane[4] - mouse.deltaPos[2]
			end
		else
			local magn = mouse.deltaPos:length() * 1000
			if magn > 0 then
				local normDelta = mouse.deltaPos / magn
				local r = quat():fromAngleAxis(-normDelta[2], normDelta[1], 0, -magn)
				if rotateClip[0] == 0 then
					viewAngle = (viewAngle * r):normalize()
				else
					local clipPlane = clipInfos[rotateClip[0]].plane
					local clipNormal = (viewAngle * r * viewAngle:conjugate()):conjugate():rotate(vec3(clipPlane:unpack()))
					for i=1,3 do
						clipPlane[i] = clipNormal[i]
					end
				end
			end
		end
	end
	
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	local w, h = self:size()
	local ar = w / h
	local znear, zfar = .1, 100
	gl.glFrustum(-ar*znear, ar*znear, -znear, znear, znear, zfar)

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	
	gl.glTranslatef(0, 0, -viewDist)

	local aa = viewAngle:toAngleAxis()
	gl.glRotated(-aa[4], aa[1], aa[2], aa[3])

	for i,clipInfo in ipairs(clipInfos) do
		gl.glClipPlane(gl.GL_CLIP_PLANE0+i-1, vec4d(clipInfo.plane:unpack()):ptr())
-- intel/ubuntu was having trouble when the clip plane included the viewport
-- so I moved the clipping code to the shader
--		if clipInfo.enabled[0] then 
--			gl.glEnable(gl.GL_CLIP_PLANE0+i-1)
--		end
	end

	gl.glTranslatef(-.5, -.5, -.5)

	volumeShader:use()
	self.solver.tex:bind(0)
	hsvTex:bind(1)
	gl.glUniform1f(volumeShader.uniforms.alpha.loc, alpha[0])
	gl.glUniform1f(volumeShader.uniforms.alphaGamma.loc, alphaGamma[0])
	gl.glUniform1iv(volumeShader.uniforms['clipEnabled[0]'].loc, 4, 
		ffi.new('int[4]', clipInfos:map(function(info) return info.enabled[0] end)))

	gl.glEnable(gl.GL_TEXTURE_GEN_S)
	gl.glEnable(gl.GL_TEXTURE_GEN_T)
	gl.glEnable(gl.GL_TEXTURE_GEN_R)
	gl.glTexGeni(gl.GL_S, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_T, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGeni(gl.GL_R, gl.GL_TEXTURE_GEN_MODE, gl.GL_OBJECT_LINEAR)
	gl.glTexGendv(gl.GL_S, gl.GL_OBJECT_PLANE, vec4d(1,0,0,0):ptr())
	gl.glTexGendv(gl.GL_T, gl.GL_OBJECT_PLANE, vec4d(0,1,0,0):ptr())
	gl.glTexGendv(gl.GL_R, gl.GL_OBJECT_PLANE, vec4d(0,0,1,0):ptr())

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
	gl.glEnable(gl.GL_BLEND)

	--[[ points
	gl.glPointSize(2)
	gl.glBegin(gl.GL_POINTS)
	for _,pt in ipairs(self.pts) do
		gl.glVertex3d( 
			(pt[1] - .5)/(self.max[1] + 1),
			(pt[2] - .5)/(self.max[2] + 1),
			(pt[3] - .5)/(self.max[3] + 1))
	end
	gl.glEnd()
	--]]
	-- [[ slices
	local n = 255
	local fwd = -viewAngle:zAxis()
	local fwddir = select(2, table(fwd):map(math.abs):sup())
	local quad = {{0,0},{1,0},{1,1},{0,1}}
	local jmin, jmax, jdir
	if fwd[fwddir] < 0 then
		jmin, jmax, jdir = 0, n, 1
	else
		jmin, jmax, jdir = n, 0, -1
	end
	gl.glUniform3f(volumeShader.uniforms.normal.loc, fwddir==1 and jdir or 0, fwddir==2 and jdir or 0, fwddir==3 and jdir or 0)
	
	gl.glBegin(gl.GL_QUADS)
	for j=jmin,jmax,jdir do
		local f = j/n
		for _,vtx in ipairs(quad) do
			if fwddir == 1 then
				gl.glVertex3f(f, vtx[1], vtx[2])
			elseif fwddir == 2 then
				gl.glVertex3f(vtx[1], f, vtx[2])
			elseif fwddir == 3 then
				gl.glVertex3f(vtx[1], vtx[2], f)
			end
		end
	end
	gl.glEnd()
	--]]
	
	gl.glDisable(gl.GL_BLEND)

	gl.glDisable(gl.GL_TEXTURE_GEN_S)
	gl.glDisable(gl.GL_TEXTURE_GEN_T)
	gl.glDisable(gl.GL_TEXTURE_GEN_R)

	hsvTex:unbind(1)
	self.solver.tex:unbind(0)
	volumeShader:useNone()

	-- TODO add back in showGradTrace / showCurl

--	for i,clipInfo in ipairs(clipInfos) do
--		gl.glDisable(gl.GL_CLIP_PLANE0+i-1)
--	end

	glreport'here'

	App.super.update(self)
end

function App:updateGUI()
	ig.igText('iteration: '..self.solver.iteration)
	
	if ig.igCombo('display', self.solver.displayVarPtr, self.solver.displayVarNames) then
		self.solver:refreshDisplayVarKernel()
	end
	ig.igText(('%.3e to %.3e'):format(self.minValue, self.maxValue))
	
	local gradImageSize = ig.ImVec2(128, 32)
	ig.igImage(
		ffi.cast('void*',ffi.cast('intptr_t',hsvTex.id)),
		gradImageSize)
	
	local gradScreenPos = ig.igGetCursorScreenPos()
	local mousePos = ig.igGetMousePos()
	local cursorX = mousePos.x - gradScreenPos.x
	local cursorY = gradScreenPos.y - mousePos.y
	if cursorX >= 0 and cursorX <= gradImageSize.x
	and cursorY >= 0 and cursorY <= gradImageSize.y
	then
		local frac = cursorX / gradImageSize.x
		ig.igBeginTooltip()
		ig.igText(tostring( self.minValue * (1-frac) + self.maxValue * frac ))
		ig.igEndTooltip()
	end

	ig.igText'transparency:'
	ig.igSliderFloat('alpha', alpha, 0, 1, '%.3e', 10)
	ig.igSliderFloat('gamma', alphaGamma, 0, 1000, '%.3e', 10)
	ig.igRadioButton("rotate camera", rotateClip, 0)
	for i,clipInfo in ipairs(clipInfos) do
		ig.igPushIDStr('clip '..i)
		ig.igCheckbox('clip', clipInfo.enabled)
		ig.igSameLine()
		ig.igRadioButton('rotate', rotateClip, i)
		ig.igPopID()
	end
	--ig.igCheckbox('show gradient trace', showGradTrace)
	--ig.igCheckbox('show curl trace', showCurlTrace)

	ig.igSeparator()
	ig.igText'simulation:'
		
	if ig.igCheckbox('converge alpha', self.solver.convergeAlpha) then
		print('alpha', self.solver.convergeAlpha[0])
		self.solver:refreshKernels()
	end
	if ig.igCheckbox('converge beta', self.solver.convergeBeta) then
		print('beta', self.solver.convergeBeta[0])
		self.solver:refreshKernels()
	end
	if ig.igCheckbox('converge gamma', self.solver.convergeGamma) then
		print('gamma', self.solver.convergeGamma[0])
		self.solver:refreshKernels()
	end

	-- InputFloat doesn't allow formats
	-- SliderFloat allows formats but doesn't allow text-editing
	-- hmm...
	local buf = ffi.new('char[256]', ('%e'):format(self.solver.updateAlpha))
	if ig.igInputText('step scale', buf, ffi.sizeof(buf)) then
		local f = tonumber(ffi.string(buf, ffi.sizeof(buf)))
		if f then
			self.solver.updateAlpha = f
		end
	end

	local bool = ffi.new('bool[1]', self.solver.useLineSearch)
	if ig.igCheckbox('line search', bool) then
		self.solver.useLineSearch = bool[0]
	end

	if ig.igButton(self.updateMethod and 'Stop' or 'Start') then
		self.updateMethod = not self.updateMethod
	end
	ig.igSameLine()
	if ig.igButton'Step' then
		self.updateMethod = 'step'
	end

	ig.igCombo('solver', self.solver.updateMethod, self.solver.updateMethods)

	ig.igSeparator()
	ig.igText'initial conditions:'
	for i,initCond in ipairs(self.solver.initConds) do
		if ig.igRadioButton(initCond.name, self.solver.initCondPtr, i) then
			self.solver:refreshInitCond()
		end
	end
	
	if ig.igButton'Reset' then
		print'resetting...'
		self.solver:resetState()
		self.updateMethod = nil
	end
end

return App
