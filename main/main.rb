require_relative '../lib/mechsimul'

include MyMath
include RGSS
RGSS.title = 'mechsimul'
rgss_main do
	Graphics.resize_screen 1024, 768
	Graphics.frame_rate = 15
	Font.default_name = 'mplus-1m-regular.ttf'
	Thread.start do
		loop do
			Graphics.update
			Input.update
		end
	end
	period_detected = false
	block = ->t, (q, p) do
		if !$canvas.detect_period && !period_detected
			period_detected = true
			Thread.start do
				IO.write 'output.txt', dft_1d(Numo::NArray.asarray($canvas.history), -1).to_a.join(?\n)
			end
		end
	end
	$canvas = Canvas.new 1,->t{t*40},->y{100*y+384},
	                     [Color.new(255,255,255),Color.new(255,255,0)],
	                     block
	$canvas.detect_period = true
	solve_hamiltonian 1,Numo::NArray[Math::PI-0.1,0.0],Float::INFINITY,
	                  1e-3,->t,(q,p){-Math.cos(q)+p**2},1e-6
end
