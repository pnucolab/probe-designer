<script>
  export let type = 'button';
  export let variant = 'primary'; // 'primary', 'secondary', 'outline', 'ghost', 'danger', 'success', 'gradient'
  export let size = 'md'; // 'sm', 'md', 'lg'
  export let disabled = false;
  export let fullWidth = false;
  export let onClick = () => {};
  
  let isHovered = false;
  let isFocused = false;

  function getButtonStyles() {
    const baseStyles = {
      border: 'none',
      fontWeight: '600',
      transition: 'all 0.2s',
      display: 'inline-flex',
      alignItems: 'center',
      justifyContent: 'center',
      gap: '8px'
    };

    const sizeStyles = {
      sm: { padding: '8px 16px', fontSize: '14px', borderRadius: '6px' },
      md: { padding: '12px 24px', fontSize: '16px', borderRadius: '8px' },
      lg: { padding: '16px 48px', fontSize: '18px', borderRadius: '12px' }
    };

    const variantStyles = {
      primary: {
        background: '#3b82f6',
        color: 'white',
        backgroundHover: '#2563eb',
        boxShadow: isHovered || isFocused ? '0 6px 12px rgba(59,130,246,0.3)' : '0 2px 4px rgba(0,0,0,0.1)'
      },
      gradient: {
        background: 'linear-gradient(135deg, #3b82f6, #2563eb)',
        color: 'white',
        backgroundHover: 'linear-gradient(135deg, #2563eb, #1d4ed8)',
        boxShadow: isHovered || isFocused ? '0 8px 16px rgba(59,130,246,0.4)' : '0 4px 8px rgba(59,130,246,0.3)',
        transform: isHovered || isFocused ? 'translateY(-3px)' : 'translateY(0)'
      },
      secondary: {
        background: '#6b7280',
        color: 'white',
        backgroundHover: '#4b5563',
        boxShadow: isHovered || isFocused ? '0 4px 8px rgba(107,114,128,0.3)' : 'none'
      },
      outline: {
        background: 'transparent',
        color: '#3b82f6',
        border: '1px solid #3b82f6',
        backgroundHover: '#3b82f6',
        colorHover: 'white'
      },
      ghost: {
        background: 'transparent',
        color: '#374151',
        backgroundHover: '#f3f4f6'
      },
      danger: {
        background: '#dc2626',
        color: 'white',
        backgroundHover: '#b91c1c',
        boxShadow: isHovered || isFocused ? '0 4px 8px rgba(220,38,38,0.3)' : 'none'
      },
      success: {
        background: '#10b981',
        color: 'white',
        backgroundHover: '#059669',
        boxShadow: isHovered || isFocused ? '0 4px 8px rgba(16,185,129,0.3)' : 'none'
      }
    };

    const variantStyle = variantStyles[variant];
    const sizeStyle = sizeStyles[size];

    let styles = {
      ...baseStyles,
      ...sizeStyle,
      background: (isHovered || isFocused) && variantStyle.backgroundHover ? variantStyle.backgroundHover : variantStyle.background,
      color: (isHovered || isFocused) && variantStyle.colorHover ? variantStyle.colorHover : variantStyle.color,
      cursor: disabled ? 'not-allowed' : 'pointer',
      opacity: disabled ? '0.5' : '1',
      width: fullWidth ? '100%' : 'auto',
      boxShadow: variantStyle.boxShadow || 'none',
      transform: variantStyle.transform || 'none',
      border: variantStyle.border || 'none'
    };

    return Object.entries(styles)
      .map(([key, value]) => `${key.replace(/[A-Z]/g, m => '-' + m.toLowerCase())}: ${value}`)
      .join('; ');
  }

  function handleClick(event) {
    if (!disabled) {
      onClick(event);
    }
  }
</script>

<button
  {type}
  {disabled}
  style={getButtonStyles()}
  on:click={handleClick}
  on:mouseover={() => isHovered = true}
  on:mouseout={() => isHovered = false}
  on:focus={() => isFocused = true}
  on:blur={() => isFocused = false}
  {...$$restProps}
>
  <slot />
</button>