require 'matrix'
require 'cmath'
require 'openrgss'
require 'bigdecimal'
require 'bigdecimal/math'

module AdditionRefinements
	[Integer, Float, BigDecimal, Complex, Rational].each { _1.prepend self }
	def + other
		zero? && [Vector, Matrix].any? { other.is_a? _1 } ? other : super(other)
	end
end

class Vector
	def self.build(...)
		new Array.new(...)
	end
	def to_ary
		to_a
	end
end

class Array
	def inner other
		zip(other).sum { _1 * _2 }
	end
end

module MyMath
	include CMath
	
	module_function
	
	def fft_recursive signal
		return signal if signal.size <= 1
		evens_odds = signal.partition.with_index { _2.even? }
		evens, odds = evens_odds.map { fft_recursive(_1) * 2 }
		coef = Complex 0, -2 * PI / signal.size
		evens.zip(odds).map.with_index do |(even, odd), i|
			even + odd * exp(coef * i)
		end
	end
	
	def stretch signal, length
		coef = signal.length / length.to_f
		Array.new length do |i|
			i *= coef
			j = i.floor
			k = j + 1
			k = 0 if k == signal.length
			signal[j] + (signal[k] - signal[j]) * (i % 1)
		end
	end
	
	def fft signal
		fft_recursive stretch signal, 2**log2(signal.size).ceil
	end
	
	FORWARD_EULER = [[],[1]]
	EXPLICIT_MIDPOINT = [[],[1/2.0],[0,1]]
	HEUN = [[],[1],[1/2.0,1/2.0]]
	RALSTON = [[],[2/3.0],[1/4.0,3/4.0]]
	KUTTA_3RD = [[],[1/2.0],[-1,2],[1/6.0,2/3.0,1/6.0]]
	HEUN_3RD = [[],[1/3.0],[0,2/3.0],[1/4.0,0,3/4.0]]
	RALSTON_3RD = [[],[1/2.0],[0,3/4.0],[2/9.0,1/3.0,4/9.0]]
	SSPRK3 = [[],[1],[1/4.0,1/4.0],[1/6.0,1/6.0,2/3.0]]
	CLASSIC_4TH = [[],[1/2.0],[0,1/2.0],[0,0,1],[1/6.0,1/3.0,1/3.0,1/6.0]]
	RALSTON_4TH = [[],[0.4],[0.29697761,0.15875964],[0.21810040,-3.05096516,3.83286476],[0.17476028, -0.55148066, 1.20553560, 0.17118478]]
	THREE_EIGHTH_4TH = [[],[1/3.0],[-1/3.0,1],[1,-1,1],[1/8.0,3/8.0,3/8.0,1/8.0]]
	
	def runge_kutta initial, max_t, dt, (*pyramid, coefs), func
		(0..max_t).step(dt).reduce initial do |ret, t|
			break unless $canvas.trace t, ret if $canvas
			coefs.zip(pyramid).each_with_object([]).sum do |(coef, row), ary|
				coef * ary.push(func.(t, row.inner(ary) * dt + ret)).last
			end * dt + ret#p(ret)
		end
	end
	
	def div x0, dx, f
		n = x0.size
		Array.new n do |i|
			dv = Vector.basis(size: n, index: i) * dx
			(f.(x0 + dv) - f.(x0 - dv)) / (2 * dx)
		end
	end
	
	def solve_hamiltonian n, qp0, max_t, dt, hamiltonian, dqp
		runge_kutta qp0, max_t, dt, CLASSIC_4TH, ->t, qp do
			dqpdt = div qp, dqp, -> { hamiltonian.(t, _1) }
			Vector[*dqpdt[n...n*2], *dqpdt[0...n].map(&:-@)]
		end
	end
end

class Canvas
	
	attr_accessor :history, :detect_period
	
	def initialize n, mapping_x, mapping_y, colors, on_trace
		@colors = colors
		@sprite = Sprite.new
		@old_sprite = Sprite.new
		@sprite.bitmap = Bitmap.new Graphics.width, Graphics.height
		@old_sprite.bitmap = Bitmap.new Graphics.width, Graphics.height
		@mapping_x = mapping_x
		@mapping_y = mapping_y
		@labels = Array.new n * 2 do |i|
			sprite = Sprite.new
			sprite.bitmap = Bitmap.new 30, 30
			sprite.bitmap.font.color = @colors[i]
			sprite.bitmap.draw_text 0, 0, 30, 30, (i < n ? "q#{i}" : "p#{i - n}"), 2
			sprite.ox = sprite.width
			sprite.oy = sprite.height / 2
			sprite.x = Graphics.width
			sprite
		end
		@quo = 0
		@detect_period = false
		@on_trace = on_trace
		@history = []
	end
	
	def trace t, ret
		if @detect_period
			if @initial
				r = (ret - @initial).r
				@last_r ||= Float::INFINITY
				if @last_r < 5e-5 && @last_r <= r
					puts "Period found."
					@detect_period = false
				end
			else
				@initial = ret
			end
			@history.push ret[0]
			@last_r = r
		end
		x = @mapping_x.(t) - Graphics.width * @quo
		@sprite.x = Graphics.width - x
		@old_sprite.x = -x
		if @sprite.x < 0
			@old_sprite, @sprite = @sprite, @old_sprite
			x -= Graphics.width
			@sprite.x = Graphics.width - x
			@sprite.bitmap.clear
			@quo += 1
		end
		ret.zip @labels, @colors do |y, label, color|
			@sprite.bitmap.set_pixel x, (label.y = @mapping_y.(y)), color
		end
		@on_trace.(t, ret)
		true
	end
	
	def visible= visibility
		@sprite.visible = visibility
		@old_sprite.visible = visibility
		@labels.each { _1.visible = visibility }
	end
end
